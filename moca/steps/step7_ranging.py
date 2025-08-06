import os
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mode as scipy_mode
from tqdm import tqdm
import re

def _custom_mode(series):
    """
    Calculates the mode of data binned into groups of 10.
    This version uses pandas' robust mode() function.
    """
    if series.empty or series.isnull().all():
        return np.nan
    
    # Floor division to bin values (e.g., 23 -> 2.0)
    binned_series = (series // 10).dropna()
    if binned_series.empty:
        return np.nan

    # pandas.Series.mode() returns a Series, as there can be multiple modes.
    # We take the first one, which matches R's which.max behavior for ties.
    mode_val = binned_series.mode()
    
    if not mode_val.empty:
        return mode_val.iloc[0] * 10
    
    return np.nan

def _calculate_stats(df, source_file, species, model_tag):
    """Groups a dataframe by motif and calculates positional statistics."""
    if df.empty: return pd.DataFrame()
    
    df = df.copy()
    df['trunc_start'] = df['start']
    
    aggregations = {
        'min': ('start', 'min'), 'max': ('start', 'max'),
        'q10': ('start', lambda x: x.quantile(0.1)),
        'median': ('start', 'median'),
        'q90': ('start', lambda x: x.quantile(0.9)),
        'mode': ('trunc_start', _custom_mode),
        'mean': ('start', 'mean'), 'sd': ('start', 'std'),
        'number': ('start', 'count')
    }
    result_df = df.groupby('motif').agg(**aggregations).reset_index()
    result_df['cv'] = (result_df['sd'] / result_df['mean']) * 100
    result_df['iqr'] = result_df['q90'] - result_df['q10']
    result_df['Species'] = species
    result_df['Model'] = model_tag
    result_df['source'] = os.path.basename(source_file)
    result_df['epm'] = "epm_" + result_df['Species'] + "_" + result_df['Model'] + "_" + result_df['motif']
    final_cols = ['epm', 'min', 'max', 'mean', 'median', 'mode', 'q10', 'q90', 'sd', 'cv', 'iqr', 'number', 'source']
    return result_df[final_cols]

def _create_motif_id(row):
    """
    Safely creates a motif ID from a dataframe row, exactly mimicking the R script's logic.
    """
    try:
        pattern_name = row['pattern']
        metacluster = row['metacluster']
        
        pattern_num_str = pattern_name[8:]
        if not pattern_num_str:
            return None
            
        pattern_num = int(pattern_num_str)
        prefix = "p0" if metacluster == "metacluster_0" else "p1"
        
        return f"{prefix}m{pattern_num:02d}"
    except (ValueError, TypeError, IndexError):
        return None

def run(config, common_settings):
    """
    Main function for positional ranging. Extracts seqlet locations from HDF5,
    calculates statistics, and generates smoothed distribution plots.
    """
    input_hdf5_path = common_settings.get('input_hdf5')
    output_dir = config['output_dir']
    if not input_hdf5_path or not os.path.exists(input_hdf5_path):
        print(f"Error: Input HDF5 file not found at '{input_hdf5_path}'. Skipping ranging analysis.")
        return

    species = common_settings.get('species_tag', 'unk')
    model = common_settings.get('model_tag', 'm0')
    range1_end = config.get('range1_end', 1500)
    range2_start = config.get('range2_start', 1520)
    sequence_length = config.get('sequence_length', 3020)
    generate_plots = config.get('generate_plots', True)

    print(f"Extracting and parsing seqlet data from {input_hdf5_path}...")
    all_seqlets_data = []
    with h5py.File(input_hdf5_path, 'r') as h5file:
        metacluster_group = h5file.get('metacluster_idx_to_submetacluster_results')
        if metacluster_group is None:
            print("Error: 'metacluster_idx_to_submetacluster_results' not found in HDF5 file.")
            return
        for mcluster_name, mcluster_data in tqdm(metacluster_group.items(), desc="Extracting Seqlets"):
            patterns_group = mcluster_data['seqlets_to_patterns_result']['patterns']
            for pattern_name, pattern_data in patterns_group.items():
                if pattern_name.startswith('pattern_'):
                    for seqlet_bytes in pattern_data['seqlets_and_alnmts']['seqlets'][:]:
                        all_seqlets_data.append({
                            'seqlets': seqlet_bytes.decode('utf-8'),
                            'pattern': pattern_name,
                            'metacluster': mcluster_name
                        })

    if not all_seqlets_data:
        print("Warning: No seqlets found in HDF5 file. Skipping ranging.")
        return

    data = pd.DataFrame(all_seqlets_data)
    pattern_regex = r'example:(\d+),start:(\d+),end:(\d+),rc:(\w+)'
    extracted_df = data['seqlets'].str.extract(pattern_regex)
    extracted_df.columns = ['example', 'start', 'end', 'rc']
    data = pd.concat([data[['pattern', 'metacluster']], extracted_df], axis=1)
    
    data['motif'] = data.apply(_create_motif_id, axis=1)
    data.dropna(subset=['motif', 'start', 'end'], inplace=True)

    data['start'] = pd.to_numeric(data['start'])
    data['end'] = pd.to_numeric(data['end'])

    if generate_plots:
        plot_dir = os.path.join(output_dir, 'distribution_plots')
        os.makedirs(plot_dir, exist_ok=True)
        print(f"Generating positional distribution plots in {plot_dir}...")
        
        for motif_id, group in data.groupby('motif'):
            plt.figure(figsize=(10, 2.5))
            sns.kdeplot(group['start'], color="darkcyan", fill=True)
            plt.title(f"EPM Positional Density for Motif: {motif_id}")
            plt.xlabel("Position")
            plt.ylabel("Density")
            plt.xlim(0, sequence_length)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            
            filename = f"epm_{species}_{model}_{motif_id}_density.png"
            plt.savefig(os.path.join(plot_dir, filename), dpi=150, bbox_inches='tight')
            plt.close()

    print("Calculating positional statistics for TSS and TTS ranges...")
    range1_df = data.query(f"start >= 1 and end <= {range1_end}").copy()
    result1 = _calculate_stats(range1_df, input_hdf5_path, species, model)
    
    range2_df = data.query(f"start >= {range2_start} and end <= {sequence_length - 20}").copy()
    range2_df['start'] = sequence_length - range2_df['start']
    result2 = _calculate_stats(range2_df, input_hdf5_path, species, model)

    if not result1.empty:
        tss_path = os.path.join(output_dir, f"{species}{model}-TSS_motif_ranges.csv")
        result1.to_csv(tss_path, index=False, float_format='%.3f')
        print(f"TSS range statistics saved to {tss_path}")

    if not result2.empty:
        tts_path = os.path.join(output_dir, f"{species}{model}-TTS_motif_ranges.csv")
        result2.to_csv(tts_path, index=False, float_format='%.3f')
        print(f"TTS range statistics saved to {tts_path}")

    print("Ranging step successfully completed.")
