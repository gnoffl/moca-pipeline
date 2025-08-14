# MOCAP: MOtif Characterization & Annotation Pipeline
MOCAP is a comprehensive and configurable bioinformatics pipeline designed to analyze DNA motifs discovered from deep learning models. It automates the entire workflow, from extracting raw motifs from modisco HDF5 files to generating a final, interactive report for key genes.

The pipeline is built in Python and designed to be modular, memory-efficient, and easily configurable through a central config.yaml file. It replaces a series of R scripts with a single, streamlined application.

## Features
**Configurable Workflow:** Control the entire pipeline from a single, human-readable .yaml file.

**Modular Steps:** Each stage of the analysis is a separate, well-defined step.

**Multiple Motif Types:** Extract and process Position Frequency Matrices (PFMs), Position Weight Matrices (PWMs), or Contribution Weight Matrices (CWMs).

**Motif Cleaning:** Optional trimming of noisy, low-information content positions from motif ends.

**Comprehensive Comparison:** Compare your motifs against multiple, curated public databases (JASPAR, Arabidopsis DAP, RNA-binding proteins, miRNAs) with separate, organized outputs.

**Positional Analysis:** Calculate detailed statistics on the genomic location of motifs relative to transcription start and end sites (TSS/TTS).

**Genomic Annotation:** Map motif occurrences to a reference genome and annotate them based on a GFF/GTF file.

**Advanced Visualization:** Generate high-quality sequence logos for all motifs and smoothed positional density plots. Produce interactive, JBrowse-like HTML reports for specific genes, complete with a detailed summary table.

**Memory Efficient:** Designed to handle large occurrence files by processing data in manageable chunks.

**Caching System:** Automatically skips steps that have already been completed to save time on re-runs.

##  Project Structure
The project is organized into a clean, standard Python package structure.

moca/
├── configs/              # All .yaml configuration files and requirements.txt
├── data/                 # Input data (HDF5 models, reference genomes, etc.)
├── results/              # All output files, organized by run name
├── scripts/              # Helper shell scripts (e.g., for blamm)
├── run_pipeline.py       # The main script to execute the entire pipeline
├── setup.py              # Makes the project installable
└── moca/                 # The core Python source code package
    ├── __init__.py
    └── steps/            # Each pipeline step is a separate script
        ├── __init__.py
        ├── step1_nomenclature.py
        └── ... (all other step scripts)

## Setup and Installation
The pipeline uses a Conda environment to manage its dependencies.

**Create the Conda Environment:**

`conda create --name moca-env python=3.9`

**Activate the Environment:**

`conda activate moca-env`

**Install Dependencies:** 
The configs/ directory should contain a requirements.txt file. Install all necessary libraries with a single command:

`pip install -r configs/requirements.txt`

**Install MOCA:** 
To make the pipeline's code accessible to Python, install it in "editable" mode from the main project directory:

`pip install -e .`

### External requirements
MOCAP was optimized and tested for the modisco output of DeepCRE (https://github.com/NAMlab/deepCRE_reimplemented/).  
MOCAP used BLAMM (https://github.com/biointec/blamm) for motif search in fasta-sequences like genomic data. Please install BLAMM and set the PATH in the config.yaml file and respective scripts-file. 
MOCAP optimizes the search of genomic sequence via partitioning. MOCAP was optimized and tested to work with samtools (https://www.htslib.org/). Please install samtools, or others (e.g. agat; https://agat.readthedocs.io/en/latest/) and set the PATH in the config.yaml file and respective scripts-file.
MOCAP uses convenient publicly available input data in addition to the modisco feature extraction. Please set up the following directory with recommended files: 
*motif_database* from https://meme-suite.org/meme/doc/download.html
*reference* this directory shall contain species specific genome (.fa, .fas) and annotation files (.gff, .gff3, .gtf)
*deepCRE_results* this directory shall contain the results of https://github.com/NAMlab/deepCRE_reimplemented, modisco results 

`cd data`
`mkdir motif_database reference deepCRE_results`
`cd ..`

## Workflow
Running the pipeline is a simple, two-step process: configure and execute.


### 1. Configure Your Analysis
All aspects of a pipeline run are controlled by a .yaml file inside the configs/ directory. You can create a new file for each analysis (e.g., configs/arabidopsis_run.yaml).

**Key Configuration Parameters:**

run_name: A unique name for your analysis. All results will be saved in a folder with this name inside results/.

input_hdf5: The path to the modisco HDF5 file you want to analyze.

species_tag & model_tag: Short identifiers for your species and model, used for naming output files.

use_trimmed_motifs: A global switch (true or false) that tells all downstream steps whether to use the raw or trimmed motif files.

Step Sections: Each step of the pipeline has its own section (e.g., nomenclature, comparison, annotation). To run a step, set run: true within its section.

### 2. Execute the Pipeline
Once your configuration file is ready, run the entire pipeline from your main moca/ directory with a single command:

`python run_pipeline.py --config configs/your_config_file.yaml`

The pipeline will execute each step marked with run: true, printing its progress to the terminal.

### Pipeline Steps Explained
Nomenclature: Extracts motifs from the HDF5 file. Can be configured to extract PFMs (for sequence patterns) or CWMs (for contribution scores). Optionally trims noisy ends from motifs.

Visualization: Generates high-quality PNG sequence logos for each motif.

Comparison: Compares your motifs against a list of local public databases and generates separate reports for each.

Importance: Calculates metadata for each motif, such as the sum of its contribution scores.

Saliency: (Optional) Visualizes saliency maps for specific genes.

Clustering: Performs hierarchical clustering of your motifs and generates a dendrogram.

Ranging: Analyzes the positional distribution of motifs and calculates statistics relative to the TSS and TTS.

Projection: Maps your motifs to a reference genome using the external blamm tool (https://github.com/biointec/blamm).

Annotation: Merges the motif occurrences with a GFF/GTF file to annotate which motifs fall within which gene regions.

Browser Viz: Generates a comprehensive, interactive HTML report for specific genes, combining a JBrowse-like gene model with a detailed summary table.

Performance: (Optional) Evaluates the predictive performance of your motifs against gene expression data.
