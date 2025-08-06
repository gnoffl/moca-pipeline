import os
import h5py
import numpy as np
import pandas as pd

def _generate_motif_name(prefix, species, model, metacluster_id, pattern_id, strand, seqlet_count):
    """Generates a standardized motif name from its components."""
    pattern_part = f"p{metacluster_id}m{pattern_id:02d}"
    return f"{prefix}_{species}_{model}_{pattern_part}{strand}_{seqlet_count}"

def run(config, common_settings):
    """
    Main function for the importance step.
    Extracts contribution score metadata (sum, max, min) for each motif
    from a modisco HDF5 file and saves it to a CSV file.
    
    Args:
        config (dict): The configuration specific to the importance step.
        common_settings (dict): The global configuration settings for the run.
    """
    # --- 1. Read configuration and set up paths ---
    input_hdf5_path = common_settings['input_hdf5']
    output_dir = config['output_dir']
    
    # Get common naming elements from global config
    prefix = common_settings.get('nomenclature_prefix', 'epm')
    species = common_settings.get('species_tag', 'unk')
    model = common_settings.get('model_tag', 'm0')
    
    # Construct output file path
    output_filename = f"{common_settings.get('date', 'nodate')}_{species}{model}_contrib_scores.csv"
    output_filepath = os.path.join(output_dir, output_filename)

    if not os.path.exists(input_hdf5_path):
        print(f"Error: Input HDF5 file not found at {input_hdf5_path}")
        return

    all_motif_stats = []

    # --- 2. Open HDF5 file and process motifs ---
    with h5py.File(input_hdf5_path, 'r') as h5file:
        metacluster_group = h5file.get('metacluster_idx_to_submetacluster_results')
        if metacluster_group is None:
            print("Error: 'metacluster_idx_to_submetacluster_results' not found in HDF5 file.")
            return

        # --- 3. Dynamically loop through metaclusters and patterns ---
        for mcluster_name, mcluster_data in metacluster_group.items():
            metacluster_id = int(mcluster_name.split('_')[-1])
            patterns_group = mcluster_data['seqlets_to_patterns_result']['patterns']
            
            for pattern_name, pattern_data in patterns_group.items():
                if not pattern_name.startswith('pattern_'):
                    continue

                pattern_id = int(pattern_name.split('_')[-1])
                
                # Get seqlet count
                seqlets = pattern_data['seqlets_and_alnmts']['seqlets']
                seqlet_count = len(seqlets)
                if seqlet_count == 0:
                    continue

                # Path to contribution scores
                contrib_scores = pattern_data.get('task0_contrib_scores')
                if contrib_scores is None:
                    continue

                # Process both forward and reverse strands
                for strand_code, strand_name in [('fwd', 'F'), ('rev', 'R')]:
                    matrix = contrib_scores[strand_code][:]
                    
                    # Generate standardized name to match other steps
                    motif_name = _generate_motif_name(
                        prefix, species, model, metacluster_id, 
                        pattern_id, strand_name, seqlet_count
                    )
                    
                    # Calculate stats
                    total_sum = np.sum(matrix)
                    max_val = np.max(matrix)
                    min_val = np.min(matrix)
                    
                    all_motif_stats.append({
                        'motif': motif_name,
                        'contrib_score_sum': total_sum,
                        'contrib_score_max': max_val,
                        'contrib_score_min': min_val,
                    })

    # --- 4. Write all stats to a single CSV file ---
    if not all_motif_stats:
        print("Warning: No motif contribution scores were extracted. Output file will not be created.")
        return
        
    print(f"Extracted stats for {len(all_motif_stats)} motifs. Saving to {output_filepath}...")
    
    results_df = pd.DataFrame(all_motif_stats)
    results_df.to_csv(output_filepath, index=False, float_format='%.6f')

    print("Importance step successfully completed.")
