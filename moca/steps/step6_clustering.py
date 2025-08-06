import os
import re
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from tqdm import tqdm

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

def _compare_pcc(pfm1, pfm2, min_overlap=5):
    """
    Calculates the max Pearson Correlation Coefficient (PCC) between two motifs.
    It uses a sliding window to find the best possible alignment.
    """
    n, m = pfm1.shape[1], pfm2.shape[1]
    max_pcc = 0.0
    
    for offset in range(-(m - min_overlap), n - min_overlap + 1):
        start1, end1 = max(0, offset), min(n, m + offset)
        start2, end2 = max(0, -offset), min(m, n - offset)
        
        if (end1 - start1) >= min_overlap:
            sub1, sub2 = pfm1[:, start1:end1], pfm2[:, start2:end2]
            flat1, flat2 = sub1.flatten(), sub2.flatten()
            if np.std(flat1) > 0 and np.std(flat2) > 0:
                pcc = np.corrcoef(flat1, flat2)[0, 1]
                if pcc > max_pcc:
                    max_pcc = pcc
    return max_pcc

def _get_newick(node, parent_dist, leaf_names, newick=''):
    """
    Converts a SciPy ClusterNode object to a Newick format string.
    """
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{parent_dist - node.dist:.4f}{newick}"
    else:
        if len(newick) > 0:
            newick = f"):{parent_dist - node.dist:.4f}{newick}"
        else:
            newick = ");"
        newick = _get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = _get_newick(node.get_right(), node.dist, leaf_names, newick=f",{newick}")
        newick = f"({newick}"
        return newick

def run(config, common_settings):
    """
    Clusters forward-strand motifs by comparing them against both orientations
    of other motifs to find the best possible similarity score.
    """
    # --- 1. Configuration ---
    use_trimmed = common_settings.get('use_trimmed_motifs', False)
    nomenclature_config = common_settings.get('nomenclature', {})
    
    if use_trimmed:
        input_filepath = nomenclature_config.get('output_file_trimmed')
        print("Clustering step will use TRIMMED motifs.")
    else:
        input_filepath = nomenclature_config.get('output_file_raw')
        print("Clustering step will use RAW motifs.")

    if not input_filepath or not os.path.exists(input_filepath):
        print(f"Error: Input JASPAR file not found. Skipping clustering.")
        return
        
    output_dir = config['output_dir']
    linkage_method = config.get('linkage_method', 'ward')

    # --- 2. Load Motifs and Separate Orientations ---
    print("Loading motifs and separating orientations...")
    all_motifs = _parse_jaspar(input_filepath)
    if not all_motifs:
        print("No motifs loaded. Aborting clustering.")
        return
    
    fwd_motifs = {mid: mdata for mid, mdata in all_motifs.items() if 'F_' in mid}
    
    # Create a map from forward to reverse complement IDs
    rc_map = {mid.replace('F_', 'R_'): mid for mid in fwd_motifs}
    rev_motifs = {rc_id: all_motifs[rc_id] for f_id, rc_id in rc_map.items() if rc_id in all_motifs}
    
    fwd_motif_ids = sorted(list(fwd_motifs.keys()))
    num_fwd_motifs = len(fwd_motif_ids)

    # --- 3. Create Comprehensive Pairwise Similarity Matrix ---
    print("Performing comprehensive pairwise comparison (Fwd vs Fwd/Rev)...")
    print("NOTE: The dissimilarity metric for clustering is (1 - Max_PCC_Score).")
    similarity_matrix = np.identity(num_fwd_motifs)

    for i in tqdm(range(num_fwd_motifs), desc="Clustering Progress"):
        for j in range(i + 1, num_fwd_motifs):
            fwd_id1 = fwd_motif_ids[i]
            fwd_id2 = fwd_motif_ids[j]
            
            pfm_f1 = fwd_motifs[fwd_id1]['matrix']
            pfm_f2 = fwd_motifs[fwd_id2]['matrix']
            
            # Find the reverse complement of the second motif
            rc_id2 = fwd_id2.replace('F_', 'R_')
            pfm_r2 = rev_motifs.get(rc_id2, {}).get('matrix')

            # Compare Fwd1 vs Fwd2
            pcc_ff = _compare_pcc(pfm_f1, pfm_f2)
            
            # Compare Fwd1 vs Rev2 (if it exists)
            pcc_fr = 0.0
            if pfm_r2 is not None:
                pcc_fr = _compare_pcc(pfm_f1, pfm_r2)
            
            # The final similarity is the best of the two comparisons
            final_pcc = max(pcc_ff, pcc_fr)
            similarity_matrix[i, j] = similarity_matrix[j, i] = final_pcc

    distance_matrix = 1 - similarity_matrix
    distance_df = pd.DataFrame(distance_matrix, index=fwd_motif_ids, columns=fwd_motif_ids)
    distance_df.to_csv(os.path.join(output_dir, "fwd_motif_distance_matrix.csv"))

    # --- 4. Perform Hierarchical Clustering on Forward Motifs ---
    print(f"Performing hierarchical clustering using '{linkage_method}' method...")
    condensed_distance = squareform(distance_matrix)
    Z = linkage(condensed_distance, method=linkage_method)

    # --- 5. Save Newick Tree ---
    tree_node = to_tree(Z, rd=False)
    tree_newick = _get_newick(tree_node, tree_node.dist, fwd_motif_ids)
    tree_filepath = os.path.join(output_dir, "fwd_motif_cluster_tree.nwk")
    with open(tree_filepath, 'w') as f:
        f.write(tree_newick)
    print(f"Newick tree for forward motifs saved to {tree_filepath}")

    # --- 6. Generate and Save Dendrogram Plot ---
    print("Generating dendrogram plot for forward motifs...")
    plt.figure(figsize=(12, max(15, num_fwd_motifs * 0.2)))
    dendrogram(
        Z,
        labels=fwd_motif_ids,
        orientation="left",
        leaf_font_size=8
    )
    plt.title(f"Forward Motif Hierarchical Clustering ({linkage_method.capitalize()} Linkage)")
    plt.xlabel("Distance (1 - Max Pearson Correlation)")
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    
    plot_path = os.path.join(output_dir, "fwd_motif_dendrogram.png")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Dendrogram plot saved to {plot_path}")

    print("Clustering step successfully completed.")
