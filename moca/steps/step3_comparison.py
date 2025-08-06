import os
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import logomaker
from tqdm import tqdm
import glob

def _parse_jaspar(filepath):
    """Parses a JASPAR-formatted file and returns a dictionary of motifs."""
    motifs = {}
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        motif_entries = re.findall(r'>.*?[\s\S]*?(?=>|$)', content)
        for entry in motif_entries:
            lines = entry.strip().split('\n')
            header = lines[0][1:].strip()
            parts = header.split(maxsplit=1)
            motif_id = parts[0]
            motif_name = parts[1] if len(parts) > 1 else motif_id
            
            matrix_rows = []
            for line in lines[1:]:
                cleaned_line = line.replace('[', '').replace(']', '').strip()
                row_parts = cleaned_line.split()
                if len(row_parts) > 1:
                    matrix_rows.append([float(x) for x in row_parts[1:]])
            
            if matrix_rows:
                matrix = np.array(matrix_rows)
                if matrix.shape[0] == 4 and matrix.ndim == 2:
                     motifs[motif_id] = {'id': motif_id, 'name': motif_name, 'matrix': matrix}
    except Exception as e:
        print(f"Error parsing JASPAR file {filepath}: {e}")
        return None
    return motifs

def _parse_meme(filepath):
    """Parses a MEME-formatted file and returns a dictionary of motifs."""
    motifs = {}
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        motif_blocks = re.split(r'\nMOTIF ', content)[1:]
        
        for block in motif_blocks:
            lines = block.strip().split('\n')
            header_parts = lines[0].split()
            motif_id = header_parts[0]
            motif_name = ' '.join(header_parts[1:]) if len(header_parts) > 1 else motif_id
            
            matrix_lines = []
            in_matrix_section = False
            for line in lines:
                if 'letter-probability matrix' in line:
                    in_matrix_section = True
                    continue
                if in_matrix_section and line.strip() and not line.startswith('URL'):
                    row_data = [float(x) for x in line.strip().split()]
                    matrix_lines.append(row_data)

            if matrix_lines:
                matrix = np.array(matrix_lines).T
                if matrix.shape[0] == 4 and matrix.ndim == 2:
                    motifs[motif_id] = {'id': motif_id, 'name': motif_name, 'matrix': matrix}
    except Exception as e:
        print(f"Error parsing MEME file {filepath}: {e}")
        return None
    return motifs

def _pfm_to_pwm(pfm, pseudocount=1.0):
    """Converts a Position Frequency Matrix to a Position Weight Matrix."""
    pfm_pseudo = np.maximum(pfm, 0) + pseudocount
    col_sums = np.sum(pfm_pseudo, axis=0, keepdims=True)
    col_sums[col_sums == 0] = 1
    return pfm_pseudo / col_sums

def _compare_motifs_pcc_sliding(motif1_pfm, motif2_pfm, min_overlap=5):
    """
    Finds the best alignment of two motifs using a sliding window and Pearson correlation.
    """
    n, m = motif1_pfm.shape[1], motif2_pfm.shape[1]
    best_result = {'pcc': -1.0, 'p_value': 1.0, 'overlap': 0, 'distance': float('inf')}
    
    if n == 0 or m == 0: return best_result

    max_pcc = -1.1

    for offset in range(-(m - min_overlap), n - min_overlap + 1):
        start1, end1 = max(0, offset), min(n, m + offset)
        start2, end2 = max(0, -offset), min(m, n - offset)
        
        overlap_len = end1 - start1
        if overlap_len >= min_overlap:
            sub1 = motif1_pfm[:, start1:end1]
            sub2 = motif2_pfm[:, start2:end2]
            
            flat1, flat2 = sub1.flatten(), sub2.flatten()
            if np.std(flat1) > 0 and np.std(flat2) > 0:
                try:
                    pcc, p_val = pearsonr(flat1, flat2)
                    if not np.isnan(pcc) and pcc > max_pcc:
                        max_pcc = pcc
                        pwm1, pwm2 = _pfm_to_pwm(sub1), _pfm_to_pwm(sub2)
                        distance = np.linalg.norm(pwm1.flatten() - pwm2.flatten())
                        best_result = {'pcc': pcc, 'p_value': p_val, 'overlap': overlap_len, 'distance': distance}
                except ValueError:
                    continue
    return best_result

def run(config, common_settings):
    """Compares motifs against a list of local public databases."""
    use_trimmed = common_settings.get('use_trimmed_motifs', False)
    nomenclature_config = common_settings.get('nomenclature', {})
    
    if use_trimmed:
        query_filepath = nomenclature_config.get('output_file_trimmed')
        print("Comparison step will use TRIMMED motifs.")
        if not query_filepath:
             print("Warning: 'use_trimmed_motifs' is true, but no trimmed file was generated. Falling back to raw motifs.")
             query_filepath = nomenclature_config.get('output_file_raw')
    else:
        query_filepath = nomenclature_config.get('output_file_raw')
        print("Comparison step will use RAW motifs.")

    if not query_filepath or not os.path.exists(query_filepath):
        print("Error: Input motif file not found. Skipping comparison.")
        return
        
    main_output_dir = config['output_dir']
    
    local_db_path = config.get('local_db_path')
    databases_to_compare = config.get('databases_to_compare', [])
    min_overlap = config.get('min_overlap', 5)
    top_n = config.get('top_n_plots', 3)
    skip_plots = config.get('skip_plots', False)
    
    if not local_db_path or not os.path.exists(local_db_path):
        print(f"Error: Local database path not found at '{local_db_path}'. Skipping comparison.")
        return

    print("Loading query motifs...")
    query_motifs = _parse_jaspar(query_filepath)
    if not query_motifs:
        print("Could not load query motifs. Aborting comparison.")
        return

    # --- Main loop to process each database individually ---
    for db_name in databases_to_compare:
        print(f"\n--- Processing database: {db_name} ---")
        db_subdir = os.path.join(local_db_path, db_name)
        if not os.path.isdir(db_subdir):
            print(f"Warning: Subdirectory '{db_name}' not found in '{local_db_path}'. Skipping.")
            continue

        db_output_dir = os.path.join(main_output_dir, db_name)
        os.makedirs(db_output_dir, exist_ok=True)
        
        # Find all .meme files in the subdirectory
        meme_files = glob.glob(os.path.join(db_subdir, '*.meme'))
        if not meme_files:
            print(f"Warning: No .meme files found in '{db_subdir}'. Skipping.")
            continue

        db_motifs = {}
        for meme_file in meme_files:
            print(f"--> Parsing {os.path.basename(meme_file)}")
            parsed = _parse_meme(meme_file)
            if parsed:
                db_motifs.update(parsed)

        if not db_motifs:
            print(f"Could not load any motifs from {db_name}. Skipping.")
            continue

        all_results = []
        for q_id, q_motif in tqdm(query_motifs.items(), desc=f"Comparing vs {db_name}"):
            for db_id, db_motif in db_motifs.items():
                result = _compare_motifs_pcc_sliding(q_motif['matrix'], db_motif['matrix'], min_overlap)
                if result['pcc'] > -1.0:
                    all_results.append({
                        'Query_ID': q_id, 'DB_ID': db_id, 'DB_Name': db_motif['name'],
                        'PCC': result['pcc'], 'P_Value': result['p_value'],
                        'Distance': result['distance'], 'Overlap': result['overlap']
                    })

        if not all_results:
            print(f"No valid comparisons were made for database {db_name}.")
            continue

        results_df = pd.DataFrame(all_results).sort_values(by=['Query_ID', 'PCC'], ascending=[True, False])
        results_df.to_csv(os.path.join(db_output_dir, 'comparison_results.tsv'), sep='\t', index=False, float_format='%.4f')
        print(f"Comparison results for {db_name} saved to {db_output_dir}")

        if not skip_plots:
            plot_dir = os.path.join(db_output_dir, 'comparison_plots')
            os.makedirs(plot_dir, exist_ok=True)
            
            for q_id, group in tqdm(results_df.groupby('Query_ID'), desc=f"Generating Plots for {db_name}"):
                top_matches = group.head(top_n)
                if top_matches.empty: continue
                
                n_plots = 1 + len(top_matches)
                fig, axes = plt.subplots(n_plots, 1, figsize=(10, 2.5 * n_plots), gridspec_kw={'hspace': 0.8})
                if n_plots == 1: axes = [axes]

                q_pwm = _pfm_to_pwm(query_motifs[q_id]['matrix'])
                q_df = pd.DataFrame(q_pwm.T, columns=['A', 'C', 'G', 'T'])
                logomaker.Logo(q_df, ax=axes[0], color_scheme='classic')
                axes[0].set_title(f"Query: {query_motifs[q_id]['name']}", loc='left')

                for i, (_, row) in enumerate(top_matches.iterrows()):
                    ax = axes[i+1]
                    db_id = row['DB_ID']
                    db_pwm = _pfm_to_pwm(db_motifs[db_id]['matrix'])
                    db_df = pd.DataFrame(db_pwm.T, columns=['A', 'C', 'G', 'T'])
                    logomaker.Logo(db_df, ax=ax, color_scheme='classic')
                    title = (f"Match {i+1}: {db_motifs[db_id]['name']} ({db_id})\n"
                             f"PCC: {row['PCC']:.3f} (p={row['P_Value']:.2e}), "
                             f"Dist: {row['Distance']:.3f}, Overlap: {int(row['Overlap']) if not pd.isna(row['Overlap']) else 0}")
                    ax.set_title(title, loc='left', fontsize=9)

                safe_filename = re.sub(r'[\\/*?:"<>|]', "_", q_id)
                plt.savefig(os.path.join(plot_dir, f"{safe_filename}_comparison.png"), dpi=150, bbox_inches='tight')
                plt.close(fig)
            print(f"Comparison plots for {db_name} saved to {plot_dir}")
    
    print("\nComparison step successfully completed for all databases.")
