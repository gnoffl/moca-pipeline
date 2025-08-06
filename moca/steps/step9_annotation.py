import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
import re

def _standardize_chr_name(series):
    """Removes 'chr' prefixes and converts to string for consistent merging."""
    return series.astype(str).str.lower().str.replace('^chr', '', regex=True)

def _load_gff(filepath):
    """Loads a GFF/GTF file, parsing it into a pandas DataFrame."""
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names, low_memory=False)
        
        df['gene_id'] = df['attributes'].str.extract(r'gene_id[= ]"([^"]+)"', expand=False)
        if df['gene_id'].isnull().all():
            df['gene_id'] = df['attributes'].str.extract(r'ID=([^;]+)', expand=False)
        
        df.dropna(subset=['gene_id'], inplace=True)
        return df
    except Exception as e:
        print(f"Error loading or parsing GFF/GTF file {filepath}: {e}")
        return None

def run(config, common_settings):
    """
    Main function for the annotation step.
    Filters motif occurrences based on genomic context and positional preferences.
    """
    # --- 1. Configuration ---
    output_dir = config['output_dir']
    
    projection_dir = common_settings.get('projection', {}).get('output_dir')
    ranging_dir = common_settings.get('ranging', {}).get('output_dir')
    ref_gff_file = config.get('reference_gff')
    
    occurrence_files = sorted(glob.glob(os.path.join(projection_dir, 'occurrences_part_*')))
    
    species_tag = common_settings.get('species_tag', 'unk')
    model_tag = common_settings.get('model_tag', 'm0')
    tss_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TSS_motif_ranges.csv")
    tts_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TTS_motif_ranges.csv")

    flank_size = common_settings.get('projection', {}).get('flank_size', 1000)
    chunk_size = config.get('chunk_size', 1000000)
    filter_method = config.get('filter_method', 'q1q9').lower()

    if not all(os.path.exists(f) for f in [tss_ranges_file, tts_ranges_file, ref_gff_file]) or not occurrence_files:
        print("Error: One or more required input files not found. Check paths and previous steps. Skipping annotation.")
        return

    # --- STAGE 1: Pre-filter all occurrence files ---
    print("Stage 1: Pre-filtering all occurrence files...")
    pre_filtered_dfs = []
    
    for i, occ_file_part in enumerate(occurrence_files):
        print(f"--> Reading file part {i+1}/{len(occurrence_files)}: {os.path.basename(occ_file_part)}")
        
        occ_iterator = pd.read_csv(
            occ_file_part, sep='\t', header=None, low_memory=False,
            names=['loc', 'source', 'motif', 'mstart', 'mend', 'score', 'strand', 'pval', 'seq'],
            chunksize=chunk_size
        )

        for occ_chunk in tqdm(occ_iterator, desc="Processing Chunks"):
            if occ_chunk.empty: continue

            loc_parts = occ_chunk['loc'].str.extract(r'([^:]+):(\d+)-(\d+)')
            occ_chunk['gene_start_flank'] = pd.to_numeric(loc_parts[1])
            occ_chunk['gene_end_flank'] = pd.to_numeric(loc_parts[2])

            filtered_chunk = occ_chunk[
                (occ_chunk['mstart'] <= flank_size) | 
                (abs(occ_chunk['mstart'] - (occ_chunk['gene_end_flank'] - occ_chunk['gene_start_flank']) + 1) <= flank_size)
            ].copy()
            
            if not filtered_chunk.empty:
                pre_filtered_dfs.append(filtered_chunk[['loc', 'motif', 'mstart', 'mend', 'score', 'strand']])

    if not pre_filtered_dfs:
        print("Warning: No motif occurrences found after pre-filtering. Annotation step will produce no output.")
        return
        
    occ_df = pd.concat(pre_filtered_dfs, ignore_index=True)
    print(f"Completed pre-filtering. Retained {len(occ_df)} occurrences for annotation.")

    # --- STAGE 2: Annotate using a direct, memory-efficient merge ---
    print("\nStage 2: Annotating pre-filtered occurrences with GFF data...")
    gene_annot_df = _load_gff(ref_gff_file)
    gene_annot_df = gene_annot_df[gene_annot_df['type'] == 'gene'][['seqid', 'start', 'end', 'strand', 'gene_id']].copy()
    gene_annot_df.rename(columns={'seqid': 'chr', 'start': 'gene_start', 'end': 'gene_end', 'strand': 'gene_strand'}, inplace=True)
    
    occ_df['loc_std'] = _standardize_chr_name(occ_df['loc'].str.split(':', n=1).str[0]) + ':' + occ_df['loc'].str.split(':', n=1).str[1]
    gene_annot_df['loc_std'] = (_standardize_chr_name(gene_annot_df['chr']) + ':' + 
                                (gene_annot_df['gene_start'] - flank_size).astype(str) + '-' + 
                                (gene_annot_df['gene_end'] + flank_size).astype(str))

    merged_df = pd.merge(occ_df, gene_annot_df, on='loc_std', how='inner')

    if merged_df.empty:
        print("--- FATAL ERROR ---")
        print("Merge between occurrence data and GFF failed. No matching locus identifiers found.")
        return

    # --- STAGE 3: Calculate Genomic Positions and Apply Final Filters ---
    print("Calculating genomic positions and applying final filters...")
    # *** FIX: Use the original 'loc' column, which is preserved during the merge on 'loc_std' ***
    merged_df['gene_start_flank'] = merged_df['loc'].str.split(':').str[1].str.split('-').str[0].astype(int)
    merged_df['gen_mstart'] = merged_df['gene_start_flank'] + merged_df['mstart'] - 1
    merged_df['gen_mend'] = merged_df['gene_start_flank'] + merged_df['mend'] - 1
    
    dist_to_tss = abs(merged_df['gen_mstart'] - merged_df['gene_start'])
    dist_to_tts = abs(merged_df['gen_mstart'] - merged_df['gene_end'])
    merged_df['dist_transc_border'] = np.minimum(dist_to_tss, dist_to_tts)
    merged_df['region'] = np.where(dist_to_tss < dist_to_tts, 'upstream', 'downstream')
    
    is_minus_strand = merged_df['gene_strand'] == '-'
    merged_df.loc[is_minus_strand, 'region'] = np.where(merged_df.loc[is_minus_strand, 'region'] == 'upstream', 'downstream', 'upstream')

    if filter_method != 'none':
        print(f"Applying positional preference filter (method: {filter_method})...")
        tss_motifs = pd.read_csv(tss_ranges_file)
        tts_motifs = pd.read_csv(tts_ranges_file)
        
        merged_df['epm'] = merged_df['motif'].str.extract(r'(epm_.+?_p\d+m\d+)', expand=False)
        
        print("\n--- DIAGNOSTICS: Comparing 'epm' identifiers for merge ---")
        print(f"First 5 'epm' identifiers from main data:\n{merged_df['epm'].unique()[:5]}")
        print(f"\nFirst 5 'epm' identifiers from TSS ranging file:\n{tss_motifs['epm'].unique()[:5]}")
        
        merged_df = pd.merge(merged_df, tss_motifs, on='epm', how='left', suffixes=('', '_tss'))
        merged_df = pd.merge(merged_df, tts_motifs, on='epm', how='left', suffixes=('', '_tts'))
        
        if merged_df['q10'].isnull().all():
            print("\n--- WARNING: Merge with ranging files failed! ---")
            print("The 'epm' identifiers did not match. This is likely because the 'model_tag' in your config")
            print("does not match the tag used when the ranging files were created.")
            print("Proceeding without positional filtering.")
            final_df = merged_df.copy()
        else:
            is_upstream = merged_df['region'] == 'upstream'
            is_downstream = merged_df['region'] == 'downstream'
            
            if filter_method == 'q1q9':
                cond_upstream = is_upstream & (merged_df['dist_transc_border'] >= merged_df['q10']) & (merged_df['dist_transc_border'] <= merged_df['q90'])
                cond_downstream = is_downstream & (merged_df['dist_transc_border'] >= merged_df['q10_tts']) & (merged_df['dist_transc_border'] <= merged_df['q90_tts'])
            elif filter_method == 'minmax':
                cond_upstream = is_upstream & (merged_df['dist_transc_border'] >= merged_df['min']) & (merged_df['dist_transc_border'] <= merged_df['max'])
                cond_downstream = is_downstream & (merged_df['dist_transc_border'] >= merged_df['min_tts']) & (merged_df['dist_transc_border'] <= merged_df['max_tts'])
            else:
                cond_upstream = pd.Series(False, index=merged_df.index)
                cond_downstream = pd.Series(False, index=merged_df.index)

            final_df = merged_df[cond_upstream | cond_downstream].copy()
    else:
        print("Skipping positional preference filter.")
        final_df = merged_df.copy()


    print(f"\nTotal occurrences after final filter: {len(final_df)}")
    print("Generating final output files...")
    
    final_df['chr'] = final_df['loc_std'].str.split(':').str[0]
    
    final_cols = ['chr', 'gene_id', 'gene_start', 'gene_end', 'gene_strand', 'motif', 'gen_mstart', 'gen_mend', 'strand', 'score', 'region', 'dist_transc_border']
    final_df_csv = final_df[[col for col in final_cols if col in final_df.columns]]
    csv_path = os.path.join(output_dir, f"annotated_motifs_{filter_method}.csv")
    final_df_csv.to_csv(csv_path, index=False, float_format='%.3f')
    print(f"Annotated CSV saved to {csv_path}")

    if not final_df.empty:
        bed_df = final_df[['chr', 'gen_mstart', 'gen_mend', 'motif', 'score', 'strand']].copy()
        bed_df.rename(columns={'chr': '#chrom', 'gen_mstart': 'chromStart', 'gen_mend': 'chromEnd', 'motif': 'name'}, inplace=True)
        bed_path = os.path.join(output_dir, f"annotated_motifs_{filter_method}.bed")
        bed_df.to_csv(bed_path, sep='\t', index=False, header=False)
        print(f"Annotated BED file saved to {bed_path}")
    
    print("Annotation step successfully completed.")
