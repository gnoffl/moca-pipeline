import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

def _parse_jaspar_file(filepath):
    """
    Parses a JASPAR-formatted file and yields individual motifs.
    """
    with open(filepath, 'r') as f:
        content = f.read()
    
    motif_entries = re.findall(r'>.*?[\s\S]*?(?=>|$)', content)
    
    for entry in motif_entries:
        lines = entry.strip().split('\n')
        title = lines[0][1:].strip()
        
        matrix_rows = []
        for line in lines[1:]:
            cleaned_line = line.replace('[', '').replace(']', '').strip()
            parts = cleaned_line.split()
            if len(parts) > 1:
                matrix_rows.append([int(float(x)) for x in parts[1:]])
            
        if matrix_rows:
            counts = np.array(matrix_rows)
            yield title, counts

def _pfm_to_pwm(pfm, pseudocount=1.0):
    """Converts a Position Frequency Matrix to a Position Weight Matrix."""
    pfm_pseudo = np.maximum(pfm, 0) + pseudocount
    col_sums = np.sum(pfm_pseudo, axis=0, keepdims=True)
    col_sums[col_sums == 0] = 1 # Avoid division by zero
    return pfm_pseudo / col_sums

def run(config, common_settings):
    """
    Generates high-quality sequence logos using logomaker.
    """
    use_trimmed = common_settings.get('use_trimmed_motifs', False)
    nomenclature_config = common_settings.get('nomenclature', {})
    
    if use_trimmed:
        input_filepath = nomenclature_config.get('output_file_trimmed')
        print("Visualization step will use TRIMMED motifs.")
        if not input_filepath:
             print("Warning: 'use_trimmed_motifs' is true, but no trimmed file was generated. Falling back to raw motifs.")
             input_filepath = nomenclature_config.get('output_file_raw')
    else:
        input_filepath = nomenclature_config.get('output_file_raw')
        print("Visualization step will use RAW motifs.")

    if not input_filepath or not os.path.exists(input_filepath):
        print(f"Error: Input motif file not found at '{input_filepath}'. Skipping visualization.")
        return

    output_dir = config['output_dir']
    print(f"Generating logos from: {input_filepath}")

    motif_count = 0
    for title, counts in _parse_jaspar_file(input_filepath):
        try:
            if counts.shape[1] == 0:
                print(f"Skipping empty motif: {title}")
                continue

            # Convert counts (PFM) to a probability matrix (PWM) for logomaker
            pwm = _pfm_to_pwm(counts)
            pwm_df = pd.DataFrame(pwm.T, columns=['A', 'C', 'G', 'T'])

            # Create an explicit figure for each plot to ensure it's not empty
            fig, ax = plt.subplots(1, 1, figsize=(counts.shape[1] * 0.5, 3))
            
            logo = logomaker.Logo(pwm_df, ax=ax, color_scheme='classic')
            logo.style_spines(visible=False)
            ax.set_ylabel("Bits")
            ax.set_title(title, fontsize=10)
            
            safe_filename = re.sub(r'[\\/*?:"<>|]', "_", title)
            output_path = os.path.join(output_dir, f'{safe_filename}.png')
            
            # Save the figure with high resolution
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig) # Close the specific figure to free memory
            motif_count += 1
        except Exception as e:
            print(f"Could not generate logo for motif '{title}': {e}")
            
    if motif_count > 0:
        print(f"Successfully generated {motif_count} logos in: {output_dir}")
    else:
        print("Warning: No valid motifs were processed.")

    print("Visualization step successfully completed.")
