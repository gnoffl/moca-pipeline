import os
import re
import numpy as np
import h5py

def _generate_motif_name(prefix, species, model, metacluster_id, pattern_id, strand, seqlet_count):
    """Generates a standardized motif name from its components."""
    pattern_part = f"p{metacluster_id}m{pattern_id:02d}"
    return f"{prefix}_{species}_{model}_{pattern_part}{strand}_{seqlet_count}"

def _format_matrix_for_jaspar(matrix):
    """Formats a numpy matrix into a JASPAR-style string."""
    rows = ['A', 'C', 'G', 'T']
    lines = []
    for i, row_char in enumerate(rows):
        row_str = ' '.join(map(str, matrix[i, :].astype(int)))
        lines.append(f"{row_char}  [ {row_str} ]")
    return "\n".join(lines)

def _trim_motif_by_ic(pfm, threshold=0.2):
    """Trims the ends of a PFM based on information content."""
    if pfm.shape[1] == 0: return pfm
    col_sums = pfm.sum(axis=0)
    col_sums[col_sums == 0] = 1
    ppm = pfm / col_sums[np.newaxis, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        ic = 2 + np.sum(np.nan_to_num(ppm * np.log2(ppm)), axis=0)
    above_threshold_indices = np.where(ic > threshold)[0]
    if len(above_threshold_indices) == 0:
        return np.array([[] for _ in range(4)])
    start_trim, end_trim = above_threshold_indices[0], above_threshold_indices[-1]
    return pfm[:, start_trim:end_trim+1]

def run(config, common_settings):
    """
    Extracts motifs (CWM, PFM, or PWM), saves the raw version, 
    and optionally saves a trimmed version.
    """
    input_hdf5_path = common_settings['input_hdf5']
    output_dir = config['output_dir']
    prefix = common_settings.get('nomenclature_prefix', 'epm')
    species = common_settings.get('species_tag', 'unk')
    model = common_settings.get('model_tag', 'm0')
    
    matrix_type = config.get('matrix_type', 'PFM').upper()
    
    trim_motifs_flag = config.get('trim_motifs', False)
    trim_threshold = config.get('trim_threshold', 0.2)

    if not os.path.exists(input_hdf5_path):
        print(f"Error: Input HDF5 file not found at {input_hdf5_path}")
        return

    raw_motifs = {}
    print(f"Extracting '{matrix_type}' motifs from {input_hdf5_path}...")
    with h5py.File(input_hdf5_path, 'r') as h5file:
        metacluster_group = h5file.get('metacluster_idx_to_submetacluster_results')
        if metacluster_group is None:
            print("Error: 'metacluster_idx_to_submetacluster_results' not found in HDF5 file.")
            return

        for mcluster_name, mcluster_data in metacluster_group.items():
            metacluster_id = int(mcluster_name.split('_')[-1])
            patterns_group = mcluster_data['seqlets_to_patterns_result']['patterns']
            
            for pattern_name, pattern_data in patterns_group.items():
                if not pattern_name.startswith('pattern_'): continue
                pattern_id = int(pattern_name.split('_')[-1])
                seqlet_count = len(pattern_data['seqlets_and_alnmts']['seqlets'])
                if seqlet_count == 0: continue

                if matrix_type == 'CWM':
                    matrix_path = 'task0_contrib_scores'
                else:
                    matrix_path = 'sequence'
                
                for strand_code, strand_name in [('fwd', 'F'), ('rev', 'R')]:
                    matrix = pattern_data[matrix_path][strand_code][:]
                    try:
                        # *** FIX: Reshape using Fortran (column-major) order to match R's behavior ***
                        matrix = matrix.flatten().reshape((4, 14), order='F')
                    except ValueError:
                        print(f"Warning: Skipping motif {pattern_name} (cannot reshape to 4x14).")
                        continue

                    if matrix_type == 'CWM':
                        abs_matrix = np.abs(matrix)
                        max_val = np.max(abs_matrix)
                        processed_matrix = np.round(seqlet_count * (abs_matrix / max_val)) if max_val > 0 else abs_matrix
                    elif matrix_type == 'PFM':
                        processed_matrix = np.round(matrix * seqlet_count)
                    else: # PWM
                        processed_matrix = np.round(matrix * 100)

                    motif_name = _generate_motif_name(
                        prefix, species, model, metacluster_id, 
                        pattern_id, strand_name, seqlet_count
                    )
                    raw_motifs[motif_name] = {'name': motif_name, 'matrix': processed_matrix}

    raw_filename = f"{prefix}_{species}_{model}_{matrix_type.lower()}-motifs_raw.jaspar"
    raw_filepath = os.path.join(output_dir, raw_filename)
    print(f"Saving {len(raw_motifs)} raw motifs to {raw_filepath}...")
    with open(raw_filepath, 'w') as f:
        for i, (motif_id, motif_data) in enumerate(raw_motifs.items()):
            f.write(f">{motif_data['name']}\n")
            f.write(_format_matrix_for_jaspar(motif_data['matrix']))
            if i < len(raw_motifs) - 1: f.write("\n")
    config['output_file_raw'] = raw_filepath

    if trim_motifs_flag:
        print(f"Trimming motifs with IC threshold > {trim_threshold} bits...")
        trimmed_motifs = {}
        for motif_id, motif_data in raw_motifs.items():
            trimmed_matrix = _trim_motif_by_ic(motif_data['matrix'].copy(), trim_threshold)
            if trimmed_matrix.shape[1] > 0:
                trimmed_motifs[motif_id] = {'name': motif_data['name'], 'matrix': trimmed_matrix}
        
        trimmed_filename = f"{prefix}_{species}_{model}_{matrix_type.lower()}-motifs_trimmed.jaspar"
        trimmed_filepath = os.path.join(output_dir, trimmed_filename)
        print(f"Saving {len(trimmed_motifs)} trimmed motifs to {trimmed_filepath}...")
        with open(trimmed_filepath, 'w') as f:
            for i, (motif_id, motif_data) in enumerate(trimmed_motifs.items()):
                f.write(f">{motif_data['name']}\n")
                f.write(_format_matrix_for_jaspar(motif_data['matrix']))
                if i < len(trimmed_motifs) - 1: f.write("\n")
        config['output_file_trimmed'] = trimmed_filepath
    
    print("Nomenclature step successfully completed.")
