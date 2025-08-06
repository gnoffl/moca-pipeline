import os
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

def run(config, common_settings):
    """
    Main function for the saliency visualization step.
    Generates importance score plots and saliency maps for a selected gene.
    
    Args:
        config (dict): The configuration specific to the saliency step.
        common_settings (dict): The global configuration settings for the run.
    """
    # --- 1. Read configuration and set up paths ---
    saliency_scores_file = config.get('saliency_scores_h5')
    metadata_file = config.get('metadata_csv')
    target_gene_id = config.get('target_gene_id')
    output_dir = config['output_dir']
    
    # Plotting-specific configs
    plot_range_start = config.get('plot_range_start', 750)
    plot_range_end = config.get('plot_range_end', 1500)
    vline_position = config.get('vline_position', 1500)

    if not all([saliency_scores_file, metadata_file, target_gene_id]):
        print("Error: Saliency step requires 'saliency_scores_h5', 'metadata_csv', and 'target_gene_id' in config.")
        return

    if not os.path.exists(saliency_scores_file):
        print(f"Error: Saliency scores file not found at {saliency_scores_file}")
        return
    if not os.path.exists(metadata_file):
        print(f"Error: Metadata CSV file not found at {metadata_file}")
        return

    # --- 2. Load and Process Data ---
    print(f"Loading data for saliency analysis of gene: {target_gene_id}")
    
    # Load metadata
    meta_df = pd.read_csv(metadata_file)
    if 'gene_id' not in meta_df.columns:
        print("Error: Metadata CSV must contain a 'gene_id' column.")
        return
        
    # Find the index corresponding to the target gene
    target_gene_info = meta_df[meta_df['gene_id'] == target_gene_id]
    if target_gene_info.empty:
        print(f"Error: Target gene '{target_gene_id}' not found in metadata file.")
        return
    gene_index = target_gene_info.index[0]

    # Load saliency scores from HDF5
    with h5py.File(saliency_scores_file, 'r') as h5file:
        saliency_scores = h5file.get('contrib_scores')
        if saliency_scores is None:
            print("Error: 'contrib_scores' not found in HDF5 file.")
            return
        
        # Extract the full saliency matrix for the target gene
        target_saliency_matrix = saliency_scores[:, :, gene_index]
        # Calculate importance scores by summing contributions across bases
        importance_scores = np.sum(target_saliency_matrix, axis=0)

    # --- 3. Generate Importance Score Line Plot ---
    print("Generating importance score line plot...")
    
    # Prepare data for plotting
    plot_data = pd.DataFrame({
        'position': range(len(importance_scores)),
        'importance': importance_scores
    }).query(f"position >= {plot_range_start} and position <= {plot_range_end}")

    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(plot_data['position'], plot_data['importance'], color='black', linewidth=0.5)
    ax.axvline(x=vline_position, color='red', linestyle='--', linewidth=1)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    
    ax.set_title(f"Contribution Scores for {target_gene_id}")
    ax.set_xlabel("Position Index")
    ax.set_ylabel("Importance Score")
    ax.set_xlim(plot_range_start, plot_range_end)
    ax.grid(False)
    
    line_plot_path = os.path.join(output_dir, f"{target_gene_id}_importance_plot.png")
    plt.savefig(line_plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Line plot saved to {line_plot_path}")

    # --- 4. Generate Saliency Map (Sequence Logo) ---
    print("Generating saliency map...")
    
    # Subset the full matrix for the specified range
    saliency_subset = target_saliency_matrix[:, plot_range_start:plot_range_end]
    
    # Create a DataFrame for logomaker
    saliency_df = pd.DataFrame(saliency_subset.T, columns=['A', 'C', 'G', 'T'])
    saliency_df.index.name = 'pos'

    try:
        # Generate the logo
        logo = logomaker.Logo(saliency_df, color_scheme='classic', ax=plt.gca())
        logo.style_spines(visible=False)
        logo.ax.set_ylabel("Contribution")
        logo.ax.set_title(f"Saliency Map for {target_gene_id}")
        
        saliency_map_path = os.path.join(output_dir, f"{target_gene_id}_saliency_map.png")
        plt.savefig(saliency_map_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saliency map saved to {saliency_map_path}")
    except Exception as e:
        print(f"Could not generate saliency map: {e}")

    print("Saliency step successfully completed.")
