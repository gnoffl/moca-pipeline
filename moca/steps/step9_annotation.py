import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
import re

def _standardize_chr_name(series):
    """
    Removes 'chr' prefixes from chromosome names and converts to a lowercase string.
    This ensures consistent chromosome formats (e.g., '1' and 'chr1' are treated as '1').
    """
    return series.astype(str).str.lower().str.replace('^chr', '', regex=True)

def _load_gff(filepath):
    """
    Loads a GFF/GTF file, parsing it into a pandas DataFrame.
    It specifically extracts the gene_id from the attributes column for annotation purposes.
    """
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names, low_memory=False)
        
        df['gene_id'] = df['attributes'].str.extract(r'gene_id[= ]"([^"]+)"', expand=False)
        if df['gene_id'].isnull().all():
            print("Could not find 'gene_id' attribute, trying 'ID=...' format.")
            df['gene_id'] = df['attributes'].str.extract(r'ID=([^;]+)', expand=False)
        
        df['gene_id'] = df['gene_id'].fillna('N/A')
        return df
    except Exception as e:
        print(f"Error loading or parsing GFF/GTF file {filepath}: {e}")
        return None

def save_filtered_results(df, output_dir, filter_name):
    """Saves a filtered DataFrame to specifically named CSV and BED files."""
    if df.empty:
        print(f"Warning: No occurrences remained for filter '{filter_name}'. No output files will be generated for this filter.")
        return

    # Handle potential column name conflict from merge.
    if 'score_x' in df.columns:
        df = df.rename(columns={'score_x': 'score'}, inplace=False)

    print(f"\nTotal occurrences after '{filter_name}' filter: {len(df)}")
    print(f"Generating output files for '{filter_name}'...")

    final_cols = ['chr', 'gene_id', 'gene_start', 'gene_end', 'gene_strand', 'motif', 'gen_mstart', 'gen_mend', 'strand', 'score', 'region', 'dist_transc_border']
    final_df_csv = df[[col for col in final_cols if col in df.columns]]

    csv_path = os.path.join(output_dir, f"annotated_motifs_{filter_name}.csv")
    final_df_csv.to_csv(csv_path, index=False, float_format='%.3f')
    print(f"Annotated CSV saved to {csv_path}")

    bed_df = df[['chr', 'gen_mstart', 'gen_mend', 'motif', 'score', 'strand']].copy()
    bed_df.rename(columns={'chr': '#chrom', 'gen_mstart': 'chromStart', 'gen_mend': 'chromEnd', 'motif': 'name'}, inplace=True)
    bed_path = os.path.join(output_dir, f"annotated_motifs_{filter_name}.bed")
    bed_df.to_csv(bed_path, sep='\t', index=False, header=True)
    print(f"Annotated BED file saved to {bed_path}")

def run(config, common_settings):
    """
    Main function for the annotation step.
    Filters motif occurrences based on genomic context and positional preferences.
    """
    # --- 1. Configuration & Path Setup ---
    output_dir = config['output_dir']
    projection_dir = common_settings.get('projection', {}).get('output_dir')
    ranging_dir = common_settings.get('ranging', {}).get('output_dir')
    ref_gff_file = config.get('reference_gff')
    chunk_size = config.get('chunk_size', 500000)

    # --- FIXED VALUES based on the original R script's logic ---
    EXTRACTION_FLANK_SIZE = 1000
    FILTER_FLANK_SIZE = 1500
    
    print(f"--- Using FIXED settings based on R script logic ---")
    print(f"Merge Key Flank Size: {EXTRACTION_FLANK_SIZE}bp (for GFF coordinate calculation)")
    print(f"Motif Filter Flank Size: {FILTER_FLANK_SIZE}bp (for pre-filtering from sequence ends)")
    print(f"----------------------------------------------------")

    occurrence_files = sorted(glob.glob(os.path.join(projection_dir, 'occurrences_part_*')))
    species_tag = common_settings.get('species_tag', 'unk')
    model_tag = common_settings.get('model_tag', 'm0')
    tss_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TSS_motif_ranges.csv")
    tts_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TTS_motif_ranges.csv")

    # --- Prepare GFF data for merging ---
    print("Loading and preparing GFF data for merging...")
    gene_annot_df = _load_gff(ref_gff_file)
    if gene_annot_df is None: return

    gene_annot_df = gene_annot_df[gene_annot_df['type'] == 'gene'].copy()
    gene_annot_df.rename(columns={'seqid': 'chr', 'start': 'gene_start', 'end': 'gene_end', 'strand': 'gene_strand'}, inplace=True)
    
    chr_std = _standardize_chr_name(gene_annot_df['chr'])
    start_flank = (gene_annot_df['gene_start'] - EXTRACTION_FLANK_SIZE).astype(str)
    end_flank = (gene_annot_df['gene_end'] + EXTRACTION_FLANK_SIZE).astype(str)
    gene_annot_df['merge_key'] = chr_std + ':' + start_flank + '-' + end_flank
    
    # --- STAGE 1 & 2: Process occurrences in chunks to conserve memory ---
    print("\nProcessing occurrence files in chunks...")
    fully_processed_chunks = []
    occ_col_names = ['loc', 'source', 'motif', 'mstart', 'mend', 'score', 'strand', 'pval', 'seq']

    for i, occ_file_part in enumerate(occurrence_files):
        print(f"--> Reading file part {i+1}/{len(occurrence_files)}: {os.path.basename(occ_file_part)}")
        try:
            for occ_chunk in tqdm(pd.read_csv(
                occ_file_part, sep='\t', header=None, low_memory=False,
                names=occ_col_names, chunksize=chunk_size
            ), desc="Processing Chunks"):
                if occ_chunk.empty: continue

                # 1. Pre-filter the chunk immediately based on position within the extracted sequence
                loc_parts = occ_chunk['loc'].str.extract(r'[^:]+:(\d+)-(\d+)')
                seq_start = pd.to_numeric(loc_parts[0])
                seq_end = pd.to_numeric(loc_parts[1])
                seq_length = seq_end - seq_start

                is_in_start_flank = occ_chunk['mstart'] <= FILTER_FLANK_SIZE
                is_in_end_flank = (seq_length - occ_chunk['mstart']) <= FILTER_FLANK_SIZE
                
                filtered_chunk = occ_chunk[is_in_start_flank | is_in_end_flank].copy()
                if filtered_chunk.empty: continue

                # 2. Annotate the small, filtered chunk
                filtered_chunk['merge_key'] = filtered_chunk['loc']
                merged_chunk = pd.merge(filtered_chunk, gene_annot_df, on='merge_key', how='inner')
                if merged_chunk.empty: continue

                # 3. Add to list of processed chunks
                fully_processed_chunks.append(merged_chunk)

        except Exception as e:
            print(f"Error processing file {occ_file_part}: {e}")

    if not fully_processed_chunks:
        print("Warning: No matching motifs found after filtering and annotation.")
        return
        
    # --- Assemble final DataFrame from processed chunks ---
    merged_df = pd.concat(fully_processed_chunks, ignore_index=True)
    print(f"\nSuccessfully processed all chunks. Total retained occurrences: {len(merged_df)}")

    # --- STAGE 3: Calculate Genomic Positions and Apply Final Filters ---
    print("\nStage 3: Calculating genomic positions and applying final filters...")
    flank_start_coord = merged_df['loc'].str.split(':').str[1].str.split('-').str[0].astype(int)
    merged_df['gen_mstart'] = flank_start_coord + merged_df['mstart'] - 1
    merged_df['gen_mend'] = flank_start_coord + merged_df['mend'] - 1
    
    is_upstream = merged_df['gen_mend'] < merged_df['gene_start']
    is_downstream = merged_df['gen_mstart'] > merged_df['gene_end']
    
    conditions = [is_upstream, is_downstream]
    choices = ['upstream', 'downstream']
    merged_df['region'] = np.select(conditions, choices, default='intragenic')

    dist_choices = [
        merged_df['gene_start'] - merged_df['gen_mend'],
        merged_df['gen_mstart'] - merged_df['gene_end']
    ]
    merged_df['dist_transc_border'] = np.select(conditions, dist_choices, default=0)

    is_minus_strand = merged_df['gene_strand'] == '-'
    upstream_on_minus = (merged_df['region'] == 'upstream') & is_minus_strand
    downstream_on_minus = (merged_df['region'] == 'downstream') & is_minus_strand
    merged_df.loc[upstream_on_minus, 'region'] = 'downstream'
    merged_df.loc[downstream_on_minus, 'region'] = 'upstream'
    
    try:
        tss_motifs = pd.read_csv(tss_ranges_file)
        tts_motifs = pd.read_csv(tts_ranges_file)
    except FileNotFoundError:
        print(f"Error: Could not read ranging files. Aborting annotation as positional filtering cannot be performed.")
        return

    merged_df['epm'] = merged_df['motif'].str.extract(r'(epm_.+?_p\d+m\d+)', expand=False)
    
    upstream_df = merged_df[merged_df['region'] == 'upstream'].copy()
    downstream_df = merged_df[merged_df['region'] == 'downstream'].copy()
    up_merged = pd.merge(upstream_df, tss_motifs, on='epm', how='left')
    down_merged = pd.merge(downstream_df, tts_motifs, on='epm', how='left')
    
    up_merged.dropna(subset=['min', 'max', 'q10', 'q90'], inplace=True)
    down_merged.dropna(subset=['min', 'max', 'q10', 'q90'], inplace=True)

    print("\n--- Applying 'minmax' filter ---")
    up_filtered_minmax = up_merged[up_merged['dist_transc_border'].between(up_merged['min'], up_merged['max'])]
    down_filtered_minmax = down_merged[down_merged['dist_transc_border'].between(down_merged['min'], down_merged['max'])]
    final_df_minmax = pd.concat([up_filtered_minmax, down_filtered_minmax], ignore_index=True)
    save_filtered_results(final_df_minmax, output_dir, "minmax")

    print("\n--- Applying 'q1q9' filter ---")
    up_filtered_q1q9 = up_merged[up_merged['dist_transc_border'].between(up_merged['q10'], up_merged['q90'])]
    down_filtered_q1q9 = down_merged[down_merged['dist_transc_border'].between(down_merged['q10'], down_merged['q90'])]
    final_df_q1q9 = pd.concat([up_filtered_q1q9, down_filtered_q1q9], ignore_index=True)
    save_filtered_results(final_df_q1q9, output_dir, "q1q9")

    print("\nAnnotation step successfully completed.")
