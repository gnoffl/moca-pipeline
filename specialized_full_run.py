import json
import os
import subprocess
from typing import Any, Dict, Optional, Tuple
import yaml
from argparse import ArgumentParser
from run_pipeline import pipeline


def parse_args() -> str:
    parser = ArgumentParser(description="Specialized Full Run of MOCA Pipeline")
    parser.add_argument('--config', "-c", type=str, required=True, help='Path to the JSON configuration file for the pipeline run.')
    args = parser.parse_args()
    return args.config


def run_command(*args: str, output: bool = False) -> Optional[str]:
    print(f"Running command: {' '.join(args)}")
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        raise RuntimeError(f"Command failed with error message: {result.stderr}")
    if output:
        return result.stdout

def get_annotation(fasta_path: str) -> Tuple[str, str]:
    final_fasta_path = fasta_path
    gff_path = ""
    output = run_command("python", "../moca_blue/moca_blue/ref_seq/create_annotation.py", "-f", fasta_path, output=True)
    if output is None:
        raise RuntimeError("really shouldnt come to this. run_command somehow didnt return output :(")
    lines = output.strip().split('\n')
    lines = [line.strip() for line in lines]
    for line in lines:
        if line.startswith("Deduplicated FASTA file saved to "):
            final_fasta_path = line.strip().split()[-1]
        if line.startswith("Annotation file created at "):
            gff_path = line.strip().split()[-1]
    if not gff_path:
        raise RuntimeError("GFF path not found in output of create_annotation.py")
    return final_fasta_path, gff_path


def update_template_params(config: Dict[str, Any], template_config: Dict[str, Any]):
    for key, value in config.items():
        if isinstance(value, dict):
            if key not in template_config:
                template_config[key] = {}
            template_config[key].update(value)
        else:
            if not isinstance(template_config.get(key, ""), type(value)):
                raise ValueError(f"Type mismatch for key '{key}': expected {type(value)}, got {type(template_config.get(key))}")
            template_config[key] = value


def run_configs(config_path: str):
    with open(config_path, 'r') as f:
        configs = json.load(f)
    
    for config in configs:
        fasta_path = config.pop("sequence_fasta")
        fasta_path, gff_path = get_annotation(fasta_path)

        with open("configs/local_test_gernot.yaml", 'r') as f:
            template_config = yaml.safe_load(f)
        update_template_params(config, template_config)
        template_config["projection"]["reference_genome_fasta"] = fasta_path
        template_config["projection"]["reference_gff"] = gff_path
        template_config["annotation"]["reference_gff"] = gff_path

        os.makedirs(os.path.join("results", template_config["run_name"]), exist_ok=False)
        run_config_path = os.path.join("results", template_config["run_name"], "run_config.yaml")
        with open(run_config_path, 'w') as f:
            yaml.dump(template_config, f)
        pipeline(run_config_path)

if __name__ == "__main__":
    config_path = parse_args()
    run_configs(config_path)