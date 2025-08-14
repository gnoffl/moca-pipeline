# ==================================================================
# File: moca/steps/step9_annotation_specialized.py (Final Version)
# ==================================================================
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob

def _load_special_gff(filepath):
    """Loads the special GFF format and extracts the core gene ID."""
    try:
        col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names)
        
        df['gene_id'] = df['seqid'].str.extract(r'(AT\dG\d+)', expand=False)
        df.dropna(subset=['gene_id'], inplace=True)
        return df
    except Exception as e:
        print(f"Error loading or parsing special GFF file {filepath}: {e}")
        return None

def run(config, common_settings):
    """
    Specialized annotation for pre-extracted sequences.
    This version preserves the unique sequence ID for downstream steps.
    """
    # --- 1. Configuration ---
    print("--- Memory-Optimized Annotation (with Sequence ID preservation) ---")
    output_dir = config['output_dir']
    
    projection_dir = common_settings.get('projection', {}).get('output_dir')
    ranging_dir = common_settings.get('ranging', {}).get('output_dir')
    ref_gff_file = config.get('reference_gff')
    
    occurrence_files = sorted(glob.glob(os.path.join(projection_dir, 'occurrences_part_*')))
    
    species_tag = common_settings.get('species_tag', 'unk')
    model_tag = common_settings.get('model_tag', 'm0')
    tss_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TSS_motif_ranges.csv")
    tts_ranges_file = os.path.join(ranging_dir, f"{species_tag}{model_tag}-TTS_motif_ranges.csv")

    filter_method = config.get('filter_method', 'minmax').lower()

    # --- Input Validation ---
    required_files = [ref_gff_file]
    if filter_method != 'none':
        required_files.extend([tss_ranges_file, tts_ranges_file])

    if not all(os.path.exists(f) for f in required_files) or not occurrence_files:
        print("Error: One or more required input files not found. Check paths and previous steps. Skipping annotation.")
        print(f"Missing: {[f for f in required_files if not os.path.exists(f)]}")
        if not occurrence_files: print("Missing occurrence files.")
        return

    # --- STAGE 1: Load Annotation Data Once ---
    print("\nStage 1: Loading all reference annotation data...")
    gene_annot_df = _load_special_gff(ref_gff_file)
    if gene_annot_df is None or gene_annot_df.empty:
        print("Fatal Error: GFF file could not be loaded or is empty. Cannot proceed.")
        return
        
    gene_annot_df = gene_annot_df[gene_annot_df['type'] == 'gene'][['gene_id', 'strand']].copy()
    gene_annot_df.rename(columns={'strand': 'gene_strand'}, inplace=True)
    gene_annot_df.drop_duplicates(subset=['gene_id'], inplace=True)

    tss_motifs, tts_motifs = None, None
    if filter_method != 'none':
        print(f"Applying positional preference filter (method: {filter_method})...")
        tss_motifs = pd.read_csv(tss_ranges_file)
        tts_motifs = pd.read_csv(tts_ranges_file)

    # --- STAGE 2: Process Occurrence Files in Chunks ---
    print(f"\nStage 2: Processing {len(occurrence_files)} occurrence file(s) in chunks...")
    
    TSS_POS = 1500
    TTS_POS = 1520
    processed_chunks = []

    for occ_file_part in tqdm(occurrence_files, desc="Annotating Chunks"):
        try:
            occ_chunk = pd.read_csv(
                occ_file_part, sep='\t', header=None, low_memory=False,
                names=['loc', 'source', 'motif', 'mstart', 'mend', 'score', 'strand', 'pval', 'seq']
            )
            # FIX: Rename 'loc' to 'sequence_id' to preserve it
            occ_chunk.rename(columns={'loc': 'sequence_id'}, inplace=True)
            occ_chunk['gene_id'] = occ_chunk['sequence_id'].str.extract(r'(AT\dG\d+)', expand=False)
            occ_chunk.dropna(subset=['gene_id'], inplace=True)

            merged_chunk = pd.merge(occ_chunk, gene_annot_df, on='gene_id', how='inner')
            if merged_chunk.empty:
                continue

            merged_chunk['region'] = np.where(merged_chunk['mstart'] <= TSS_POS, 'upstream', 'downstream')
            merged_chunk['dist_transc_border'] = np.where(
                merged_chunk['region'] == 'upstream',
                TSS_POS - merged_chunk['mstart'],
                merged_chunk['mstart'] - TTS_POS
            )
            
            is_minus_strand = merged_chunk['gene_strand'] == '-'
            merged_chunk.loc[is_minus_strand, 'region'] = np.where(merged_chunk.loc[is_minus_strand, 'region'] == 'upstream', 'downstream', 'upstream')
            
            if filter_method != 'none':
                merged_chunk['epm'] = merged_chunk['motif'].str.extract(r'(epm_.+?_p\d+m\d+)', expand=False)
                
                up_chunk = merged_chunk[merged_chunk['region'] == 'upstream']
                down_chunk = merged_chunk[merged_chunk['region'] == 'downstream']

                up_merged = pd.merge(up_chunk, tss_motifs, on='epm', how='left')
                down_merged = pd.merge(down_chunk, tts_motifs, on='epm', how='left')

                if filter_method == 'minmax':
                    up_filtered = up_merged[up_merged['dist_transc_border'].between(up_merged['min'], up_merged['max'])]
                    down_filtered = down_merged[down_merged['dist_transc_border'].between(down_merged['min'], down_merged['max'])]
                else: # q1q9
                    up_filtered = up_merged[up_merged['dist_transc_border'].between(up_merged['q10'], up_merged['q90'])]
                    down_filtered = down_merged[down_merged['dist_transc_border'].between(down_merged['q10'], down_merged['q90'])]
                
                processed_chunk = pd.concat([up_filtered, down_filtered], ignore_index=True)
            else:
                processed_chunk = merged_chunk

            processed_chunks.append(processed_chunk)

        except Exception as e:
            print(f"Warning: Could not process chunk {os.path.basename(occ_file_part)}: {e}")
            
    if not processed_chunks:
        print("\nError: No motif occurrences remained after processing and filtering. Cannot generate output.")
        return

    # --- STAGE 3: Final Consolidation and Output ---
    print("\nStage 3: Consolidating results and saving output...")
    final_df = pd.concat(processed_chunks, ignore_index=True)
    
    print(f"Total occurrences after final filter: {len(final_df)}")
    
    # FIX: Add 'sequence_id' to the final output file
    final_cols = ['sequence_id', 'gene_id', 'gene_strand', 'motif', 'mstart', 'mend', 'strand', 'score', 'region', 'dist_transc_border']
    final_df_csv = final_df[[col for col in final_cols if col in final_df.columns]]
    csv_path = os.path.join(output_dir, f"annotated_motifs_{filter_method}.csv")
    final_df_csv.to_csv(csv_path, index=False, float_format='%.3f')
    print(f"Annotated CSV saved to {csv_path}")
    
    print("\nAnnotation step successfully completed.")
