import os
import subprocess
import shutil
import re
from tqdm import tqdm

def _run_command(command, cwd=None):
    """Executes a simple shell command and checks for errors."""
    print(f"Running command: {' '.join(command)}")
    try:
        result = subprocess.run(
            command, cwd=cwd, check=True, capture_output=True, text=True,
        )
        print("--- Stdout ---")
        print(result.stdout)
        if result.stderr:
            print("--- Stderr ---")
            print(result.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"--- FATAL ERROR during command execution ---")
        print(f"Command failed: {' '.join(e.cmd)}")
        print(f"Return code: {e.returncode}")
        print("--- Stdout ---")
        print(e.stdout)
        print("--- Stderr ---")
        print(e.stderr)
        return False

def _run_blamm_scan_with_progress(command, cwd=None):
    """Executes the blamm scan command while showing a tqdm progress bar."""
    print(f"Running command: {' '.join(command)}")
    
    process = subprocess.Popen(
        command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    pbar = tqdm(total=100, desc="BLAMM Scan Progress", unit="%")
    
    for line in iter(process.stdout.readline, ''):
        match = re.search(r"Progress...\s*([\d\.]+)%", line)
        if match:
            progress = float(match.group(1))
            pbar.n = round(progress, 1)
            pbar.refresh()
    
    pbar.close()
    
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"--- FATAL ERROR during blamm scan ---")
        print(f"Return code: {process.returncode}")
        print("--- Stdout ---")
        print(stdout)
        print("--- Stderr ---")
        print(stderr)
        return False
        
    return True


def run(config, common_settings):
    """
    Main function for the projection step.
    Prepares a reference genome, maps motifs, and splits the large output file.
    """
    output_dir = config['output_dir']
    
    samtools_path = config.get('samtools_executable', 'samtools')
    blamm_path = config.get('blamm_executable')
    ref_genome_fasta = config.get('reference_genome_fasta')
    ref_genome_gff = config.get('reference_gff')
    
    flank_size = config.get('flank_size', 1000)
    p_value_threshold = config.get('p_value_threshold', 0.001)
    
    split_occurrences = config.get('split_occurrences', True)
    split_line_count = config.get('split_line_count', 10000000)

    use_trimmed = common_settings.get('use_trimmed_motifs', False)
    nomenclature_config = common_settings.get('nomenclature', {})
    input_motifs_file = (nomenclature_config.get('output_file_trimmed') if use_trimmed 
                         else nomenclature_config.get('output_file_raw'))

    if not all([blamm_path, ref_genome_fasta, ref_genome_gff, input_motifs_file]):
        print("Error: Projection step requires 'blamm_executable', 'reference_genome_fasta', 'reference_gff', and a valid motif file.")
        return

    # --- 1. Prepare Reference Sequence ---
    print("\n--- Preparing reference sequence with flanking regions ---")
    prepared_fasta_name = f"{os.path.basename(ref_genome_fasta).split('.')[0]}_{flank_size}bp-flank.fa"
    prepared_fasta_path = os.path.join(output_dir, prepared_fasta_name)
    ranges_file = os.path.join(output_dir, "extracted_ranges.txt")

    print("Extracting gene ranges from GFF...")
    with open(ref_genome_gff, 'r') as gff, open(ranges_file, 'w') as out:
        for line in gff:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) > 8 and 'gene' in parts[2]:
                chrom, start, end = parts[0], int(parts[3]), int(parts[4])
                new_start = max(1, start - flank_size)
                new_end = end + flank_size
                out.write(f"{chrom}:{new_start}-{new_end}\n")

    print("Using samtools to extract FASTA sequences...")
    samtools_command = [samtools_path, 'faidx', ref_genome_fasta, '-r', ranges_file, '-o', prepared_fasta_path]
    if not _run_command(samtools_command):
        print("Failed to prepare reference sequence. Aborting projection.")
        return

    # --- 2. Run BLAMM Mapping ---
    print("\n--- Running BLAMM for motif mapping ---")
    blamm_dir = os.path.dirname(blamm_path)
    blamm_exec = os.path.basename(blamm_path)
    
    manifest_file_path = os.path.join(blamm_dir, 'sequences.mf')
    print(f"Creating BLAMM manifest file at: {manifest_file_path}")
    with open(manifest_file_path, 'w') as mf:
        absolute_fasta_path = os.path.abspath(prepared_fasta_path)
        mf.write(f"group1 {absolute_fasta_path}\n")

    dict_command = [f'./{blamm_exec}', 'dict', 'sequences.mf']
    if not _run_command(dict_command, cwd=blamm_dir):
        print("BLAMM dictionary creation failed.")
        return

    hist_command = [f'./{blamm_exec}', 'hist', '-e', os.path.abspath(input_motifs_file), 'sequences.mf']
    if not _run_command(hist_command, cwd=blamm_dir):
        print("BLAMM histogram creation failed.")
        return

    scan_command = [f'./{blamm_exec}', 'scan', '-rc', '-pt', str(p_value_threshold), os.path.abspath(input_motifs_file), 'sequences.mf']
    if not _run_blamm_scan_with_progress(scan_command, cwd=blamm_dir):
        print("BLAMM scan failed.")
        return

    # --- 3. Split the Large Occurrences File ---
    if split_occurrences:
        print("\n--- Splitting large occurrences file ---")
        occurrences_file_path = os.path.join(blamm_dir, 'occurrences.txt')
        if os.path.exists(occurrences_file_path):
            split_command = [
                'split', '-l', str(split_line_count), '--numeric-suffixes',
                'occurrences.txt', 'occurrences_part_'
            ]
            if _run_command(split_command, cwd=blamm_dir):
                print("Successfully split occurrences.txt into smaller parts.")
                os.remove(occurrences_file_path)
        else:
            print("Warning: occurrences.txt not found, skipping split.")

    # --- 4. Organize BLAMM Output ---
    print("\n--- Organizing BLAMM output files ---")
    blamm_output_files = ['hist_empirical.txt', 'hist_theoretical.txt', 'PWMthresholds.txt']
    for filename in blamm_output_files:
        src = os.path.join(blamm_dir, filename)
        dst = os.path.join(output_dir, filename)
        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"Moved {filename} to results directory.")
            
    for filename in os.listdir(blamm_dir):
        if filename.startswith('occurrences_part_'):
            src = os.path.join(blamm_dir, filename)
            dst = os.path.join(output_dir, filename)
            shutil.move(src, dst)
    print("Moved split occurrence files to results directory.")

    print("Projection step successfully completed.")
