/moca_blue_project/
|
|-- moca_blue/                # The core Python source code package
|   |-- __init__.py
|   |-- steps/                # Each major pipeline step is a submodule
|   |   |-- __init__.py
|   |   |-- 1_nomenclature.py # Extracts motifs from HDF5
|   |   |-- 2_clustering.py   # Clusters motifs
|   |   |-- 3_projection.py   # Prepares for mapping
|   |   |-- 4_annotation.py   # Annotates results with gene info
|   |-- utils/                # Utility functions (e.g., file I/O, HDF5 readers)
|   |   |-- __init__.py
|   |   |-- readers.py
|   |   |-- plotting.py
|
|-- configs/                  # All analysis configurations live here
|   |-- default_config.yaml
|   |-- analysis_zea_mays.yaml
|   |-- analysis_athal.yaml
|
|-- data/                     # Raw input data
|   |-- hdf5_models/
|   |-- reference_genomes/
|
|-- results/                  # All output is neatly organized by run name
|   |-- analysis_zea_mays/
|   |   |-- 1_nomenclature/
|   |   |-- 2_clustering/
|   |   |-- 3_projection/
|   |   |-- logs/
|
|-- scripts/                  # Wrapper scripts for external tools
|   |-- run_blamm.sh
|
|-- run_pipeline.py           # The SINGLE script to execute the entire pipeline
|-- requirements.txt
|-- README.md
