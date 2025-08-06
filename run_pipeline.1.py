import argparse
import yaml
import os
import sys
from datetime import datetime
import importlib

def load_config(config_path):
    """Loads and validates the YAML configuration file."""
    print(f"Loading configuration from: {config_path}")
    if not os.path.exists(config_path):
        print(f"Error: Configuration file not found at {config_path}")
        sys.exit(1)
        
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    if 'run_name' not in config:
        print("Error: Config file must contain 'run_name'.")
        sys.exit(1)
        
    return config

def setup_directories(config):
    """Creates the necessary output directories for the run."""
    run_name = config['run_name']
    results_dir = os.path.join('results', run_name)
    os.makedirs(results_dir, exist_ok=True)
    
    # The 'trimming' step is now part of 'nomenclature', so it's removed from this list
    all_steps = ['nomenclature', 'visualization', 'comparison', 'importance', 'saliency', 
                 'clustering', 'ranging', 'projection', 'annotation', 'browser_viz', 'performance', 'logs']
    
    for step_name in all_steps:
        step_dir = os.path.join(results_dir, step_name)
        os.makedirs(step_dir, exist_ok=True)
        if step_name not in config:
            config[step_name] = {}
        config[step_name]['output_dir'] = step_dir
        
    print(f"Results will be saved in: {results_dir}")
    return config

def run_step(step_name, step_number, total_steps, config, common_settings):
    """Dynamically imports and runs a pipeline step, with caching."""
    if not config.get(step_name, {}).get('run'):
        print(f"\n[Step {step_number}/{total_steps}] Skipping {step_name.capitalize()} (disabled in config).")
        return

    # --- Caching Logic ---
    force_rerun = common_settings.get('force_rerun', False)
    output_dir = config[step_name]['output_dir']
    
    expected_outputs = {
        'nomenclature': [os.path.join(output_dir, f"{common_settings.get('nomenclature_prefix', 'epm')}_{common_settings.get('species_tag', 'unk')}_{common_settings.get('model_tag', 'm0')}_{config.get('nomenclature', {}).get('matrix_type', 'PFM').lower()}-motifs_raw.jaspar")],
        'ranging': [os.path.join(output_dir, f"{common_settings.get('species_tag', 'unk')}{common_settings.get('model_tag', 'm0')}-TSS_motif_ranges.csv")],
        'projection': [os.path.join(output_dir, 'occurrences_part_00')],
        'annotation': [os.path.join(output_dir, f"annotated_motifs_{config.get('annotation', {}).get('filter_method', 'q1q9')}.csv")]
    }

    if not force_rerun and step_name in expected_outputs:
        if step_name == 'nomenclature' and config.get('nomenclature', {}).get('trim_motifs'):
            prefix, species, model, matrix_type = (common_settings.get(k, v) for k, v in [('nomenclature_prefix', 'epm'), ('species_tag', 'unk'), ('model_tag', 'm0'), ('matrix_type', 'PFM')])
            expected_outputs['nomenclature'].append(os.path.join(output_dir, f"{prefix}_{species}_{model}_{matrix_type.lower()}-motifs_trimmed.jaspar"))

        if all(os.path.exists(f) for f in expected_outputs[step_name]):
            print(f"\n[Step {step_number}/{total_steps}] Skipping {step_name.capitalize()} (output files already exist).")
            if step_name == 'nomenclature':
                config['nomenclature']['output_file_raw'] = expected_outputs['nomenclature'][0]
                if config.get('nomenclature', {}).get('trim_motifs'):
                    config['nomenclature']['output_file_trimmed'] = expected_outputs['nomenclature'][1]
            return

    # --- Dependency Check ---
    dependencies = {
        'visualization': 'nomenclature',
        'comparison': 'nomenclature',
        'clustering': 'nomenclature',
        'projection': 'nomenclature',
        'annotation': ['projection', 'ranging'],
        'browser_viz': 'annotation',
        'performance': 'annotation'
    }
    
    deps_to_check = dependencies.get(step_name, [])
    if isinstance(deps_to_check, str): deps_to_check = [deps_to_check]

    for dep in deps_to_check:
        if not common_settings.get(dep, {}).get('run'):
            print(f"\n[Step {step_number}/{total_steps}] Skipping {step_name.capitalize()}: Dependency '{dep}' was not run.")
            return

    print(f"\n[Step {step_number}/{total_steps}] Running {step_name.capitalize()}...")
    
    try:
        module_name = f"moca.steps.step{step_number}_{step_name}"
        step_module = importlib.import_module(module_name)
        step_module.run(config=config[step_name], common_settings=common_settings)

    except ImportError as e:
        print(f"--- FATAL ERROR ---")
        print(f"Could not import the module for step: '{step_name}'.")
        print(f"Please ensure the file 'moca/steps/step{step_number}_{step_name}.py' exists and is correct.")
        print(f"Details: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"--- FATAL ERROR ---")
        print(f"An error occurred while running the '{step_name}' step.")
        print(f"Details: {e}")
        sys.exit(1)

def main():
    """Main function to orchestrate the MOCA pipeline."""
    parser = argparse.ArgumentParser(description="MOCA: Motif Characterization & Annotation Pipeline")
    parser.add_argument('--config', type=str, required=True, help='Path to the YAML configuration file for the pipeline run.')
    args = parser.parse_args()

    config = load_config(args.config)
    if 'date' not in config:
        config['date'] = datetime.now().strftime('%Y%m%d')
    config = setup_directories(config)

    print(f"\n--- Starting Pipeline Run: {config['run_name']} ---")
    start_time = datetime.now()

    # Updated list of steps in the correct order
    pipeline_steps = [
        'nomenclature', 'visualization', 'comparison', 'importance', 
        'saliency', 'clustering', 'ranging', 'projection', 
        'annotation', 'browser_viz', 'performance'
    ]
    
    total_steps = len(pipeline_steps)
    for i, step_name in enumerate(pipeline_steps):
        run_step(step_name, i + 1, total_steps, config, config)

    end_time = datetime.now()
    print(f"\n--- Pipeline Run Finished ---")
    print(f"Total execution time: {end_time - start_time}")

if __name__ == "__main__":
    main()
