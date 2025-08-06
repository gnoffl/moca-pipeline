import argparse
import yaml
import os
import sys
from datetime import datetime
import importlib.util

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
    
    all_steps = ['nomenclature', 'visualization', 'comparison', 'importance', 'saliency', 
                 'clustering', 'ranging', 'projection', 'annotation', 'performance', 'logs']
    
    for step_name in all_steps:
        step_dir = os.path.join(results_dir, step_name)
        os.makedirs(step_dir, exist_ok=True)
        if step_name not in config:
            config[step_name] = {}
        config[step_name]['output_dir'] = step_dir
        
    print(f"Results will be saved in: {results_dir}")
    return config

def run_step(step_name, step_number, total_steps, config, common_settings):
    """Dynamically imports and runs a pipeline step."""
    if not config.get(step_name, {}).get('run'):
        print(f"\n[Step {step_number}/{total_steps}] Skipping {step_name.capitalize()}.")
        return

    # Check for dependencies
    dependencies = {
        'visualization': 'nomenclature',
        'comparison': 'nomenclature',
        'clustering': 'nomenclature',
        'projection': 'nomenclature',
        'annotation': ['projection', 'ranging'],
        'performance': 'annotation'
    }
    
    deps_to_check = dependencies.get(step_name, [])
    if isinstance(deps_to_check, str):
        deps_to_check = [deps_to_check]

    for dep in deps_to_check:
        if not common_settings.get(dep, {}).get('run'):
            print(f"\n[Step {step_number}/{total_steps}] Skipping {step_name.capitalize()}: Dependency '{dep}' was not run.")
            return

    print(f"\n[Step {step_number}/{total_steps}] Running {step_name.capitalize()}...")
    
    try:
        # Construct the module name and path
        module_name = f"moca.steps.step{step_number}_{step_name}"
        
        # Dynamically import the module
        step_module = importlib.import_module(module_name.replace(f"step{step_number}", f"step{step_number-1}", 1) if step_number > 9 else module_name)
        
        # Run the 'run' function from the imported module
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
    parser.add_argument(
        '--config', 
        type=str, 
        required=True, 
        help='Path to the YAML configuration file for the pipeline run.'
    )
    args = parser.parse_args()

    config = load_config(args.config)
    if 'date' not in config:
        config['date'] = datetime.now().strftime('%Y%m%d')
    config = setup_directories(config)

    print(f"\n--- Starting Pipeline Run: {config['run_name']} ---")
    start_time = datetime.now()

    # Define the pipeline steps in order
    pipeline_steps = [
        'nomenclature', 'visualization', 'comparison', 'importance', 
        'saliency', 'clustering', 'ranging', 'projection', 
        'annotation', 'performance'
    ]
    
    total_steps = len(pipeline_steps)
    for i, step_name in enumerate(pipeline_steps):
        run_step(step_name, i + 1, total_steps, config, config)

    end_time = datetime.now()
    print(f"\n--- Pipeline Run Finished ---")
    print(f"Total execution time: {end_time - start_time}")


if __name__ == "__main__":
    main()
