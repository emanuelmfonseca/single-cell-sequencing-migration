# Load the config file
configfile: "config.yaml"

# Rule to handle the full workflow, expanding to cover all runs
rule all:
    input:
        #expand("data/run_{run_id}/fastqs", run_id=range(1, config["weekly_runs"] + 1)),
        #reference = "data/human_reference",
        normalized_h5ad="data/merged_runs/merged_adata_normalized.h5ad",
        annotate_h5ad = "data/merged_runs/merged_adata_annotate_h5ad.h5ad",
        qc_plot_unfiltered = "data/merged_runs/figures/qc_metrics_plot_unfiltered.png",
        qc_plot_filtered = "data/merged_runs/figures/qc_metrics_plot_filtered.png",
        umap = "data/merged_runs/figures/umap_merged.png",
        t_sne_plot = "data/merged_runs/figures/t_sne_merged.png",
        clustering = "data/merged_runs/clustering",
        de = "data/merged_runs/figures/rank_genes_groups_cell_type_merged.png",

# Rule to download raw sequencing data for a specific run
rule get_raw_data:
    output:
        # Define the directory where the extracted FASTQ files will be stored
        fastqs = directory("data/run_{run_id}/fastqs")
    params:
        # URL to the dataset to be downloaded, defined in the configuration file
        url = config["sc_dataset_url"],
        # Path where the downloaded tar file will be saved before extraction
        tar_file = "data/run_{run_id}/" + config["sc_dataset"]
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        """
        # Create the output directory if it doesn't exist
        mkdir -p {output.fastqs}
        
        # Download the FASTQ dataset tarball from the specified URL
        wget -N -O {params.tar_file} {params.url}
        
        # Extract the contents of the tarball into the output directory
        tar --strip-components=1 -xvf {params.tar_file} -C {output.fastqs}
                
        # Remove the tarball after extraction to save disk space
        rm {params.tar_file}
        """


# Rule to download and extract the human reference dataset
rule get_reference:
    output:
        # Define the directory where the reference files will be stored
        reference = directory("data/human_reference")
    params:
        # URL to the reference dataset tarball, defined in the configuration file
        url = config["ref_dataset_url"],
        # Path where the tar file will be saved before extraction
        tar_file = "data/human_reference/" + config["ref_dataset"]
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        """
        # Create the output directory if it doesn't already exist
        mkdir -p {output.reference}
        
        # Download the reference tarball from the specified URL
        wget -N -O {params.tar_file} {params.url}
        
        # Extract the tarball into the specified directory
        tar --strip-components=1 -xvf {params.tar_file} -C {output.reference}     
        
        # Remove the tarball after extraction to save disk space
        rm {params.tar_file}
        """


# Rule to run Cell Ranger on a set of FASTQ files
rule run_cell_ranger:
    input:
        # Directory containing the human reference genome for Cell Ranger
        reference_dir = "data/human_reference/",
        # Directory containing the FASTQ files, with run_id being a dynamic parameter
        fastqs_dir = "data/run_{run_id}/fastqs/",
    output:
        # Output directory for Cell Ranger results; it expects to generate multiple files within this directory
        # Use 'directory()' to indicate that this is an entire directory, not just a file
        cellranger = directory("data/run_{run_id}/cellranger/outs/"),  
    params:
        # Parameters extracted from the configuration file: weekly_runs (to track or group runs)
        weekly_run = config["weekly_runs"],
        # Expected number of cells, which will be used by Cell Ranger to optimize its cell-counting algorithm
        expected_cells = config["expected_cells"],
        # Number of cores to be used for the local computation
        ncores = config["ncores"],
        # Path to the installed version of Cell Ranger
        cell_ranger = config["cell_ranger_path"],
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Shell command to run Cell Ranger 'count' function.
        """
        # Extract sample prefix by finding FASTQ file names and removing _S1.* part
        sample_prefix=$(ls {input.fastqs_dir}/.fastq.gz | sed -n 's/^\(.*\)_S1.*/\1/p' | uniq)
        
        # Run the Cell Ranger 'count' tool using the extracted sample prefix, reference genome, and other params
        {params.cell_ranger}/cellranger count \
        --id cellranger \  # Identifier for the run
        --transcriptome {input.reference_dir} \  # Reference genome directory
        --fastqs {input.fastqs_dir} \  # Input FASTQ files directory
        --sample $sample_prefix \  # Sample prefix extracted from FASTQ filenames
        --expected-cells {params.expected_cells} \  # Number of expected cells
        --localcores {params.ncores}  # Number of cores for local computation
        """


# Rule to load and process data using Scanpy from the 10X Cell Ranger matrix
rule load_data:
    input:
        # Path to the script that processes the Cell Ranger matrix
        script="scripts/load_data.py",
        # Path to the Cell Ranger output matrix (filtered feature barcode matrix)
        cellranger_matrix = "data/run_{run_id}/cellranger/outs/filtered_feature_bc_matrix"
    output:
        # Define the path where the processed AnnData (H5AD) file will be saved
        h5ad = "data/merged_runs/individual_runs/run{run_id}_adata.h5ad"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        "python {input.script} {input.cellranger_matrix} {output}"


# Rule to merge individual AnnData (H5AD) files into a single merged file
rule merge_runs:
    input:
        # Path to the Python script responsible for merging the runs
        script="scripts/merge_runs.py",
        # Expanded list of individual AnnData (H5AD) files from different runs
        h5ad = expand("data/merged_runs/individual_runs/run{run_id}_adata.h5ad", run_id=range(1, config["weekly_runs"] + 1)),
    output:
        # Output path for the final merged AnnData file
        merged_h5ad = "data/merged_runs/merged_data.h5ad"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {output} {input.h5ad}"


# Rule to plot quality control metrics like percent.mt, nFeature_RNA, and nCount_RNA
rule qc_metrics_plot_unfiltered:
    input:
        # Path to the Python script that generates quality control metrics
        script="scripts/metrics_cell_qc.py",
        # Input file: the filtered AnnData file
        merged_h5ad = "data/merged_runs/merged_data.h5ad"
    output:
        # Output file: the QC plot (for example, as PNG)
        qc_plot = "data/merged_runs/figures/qc_metrics_plot_unfiltered.png"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {input.merged_h5ad} {output}"

# Rule to filter the merged AnnData (H5AD) file
rule filter_data:
    input:
        # Path to the Python script responsible for filtering the data
        script="scripts/filter_data.py",
        # Input file: the merged AnnData file
        merged_h5ad = "data/merged_runs/merged_data.h5ad"
    output:
        # Output path for the filtered AnnData file
        filtered_h5ad = "data/merged_runs/merged_adata_filtered.h5ad"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {input.merged_h5ad} {output}"

# Rule to plot quality control metrics like percent.mt, nFeature_RNA, and nCount_RNA
rule qc_metrics_plot_filtered:
    input:
        # Path to the Python script that generates quality control metrics
        script="scripts/metrics_cell_qc.py",
        # Input file: the filtered AnnData file
        filtered_h5ad = "data/merged_runs/merged_adata_filtered.h5ad"
    output:
        # Output file: the QC plot (for example, as PNG)
        qc_plot = "data/merged_runs/figures/qc_metrics_plot_filtered.png"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {input.filtered_h5ad} {output}"


# Rule to normalize the filtered AnnData (H5AD) file
rule normalize_data:
    input:
        # Path to the Python script responsible for normalizing the data
        script="scripts/normalize_data.py",
        # Input file: the hvg AnnData file
        filtered_h5ad = "data/merged_runs/merged_adata_filtered.h5ad"
    output:
        # Output path for the normalized AnnData file
        normalized_h5ad = "data/merged_runs/merged_adata_normalized.h5ad"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {input.filtered_h5ad} {output}"


# Rule to annotate the filtered AnnData (H5AD) file
rule annotate_data:
    input:
        # Python script for annotating the data
        script="scripts/annotate_cells.py",
        # Input: filtered AnnData file after high-variable gene (HVG) selection
        filtered_h5ad = "data/merged_runs/merged_adata_filtered.h5ad"
    output:
        # Output: annotated AnnData file
        annotate_h5ad = "data/merged_runs/merged_PBMC_cell_type_annotations.csv"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Run the annotation script with the filtered AnnData file as input
        "python {input.script} {input.filtered_h5ad} {output.annotate_h5ad}"


# Rule to filter the merged AnnData (H5AD) file
rule filter_hvg:
    input:
        # Path to the Python script responsible for keeping highly variable genes
        script="scripts/filter_hvg.py",
        # Input file: the merged AnnData file
        normalized_h5ad = "data/merged_runs/merged_adata_normalized.h5ad"
    output:
        # Output path for the filtered AnnData file and plot
        hvg_h5ad = "data/merged_runs/merged_adata_hvg.h5ad",
        hvg_plot = "data/merged_runs/figures/hvg_merged.png",

    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {input.normalized_h5ad} {output.hvg_h5ad} {output.hvg_plot}"


# Rule to normalize the highly variable gene (HVG) filtered AnnData (H5AD) file
rule scaling_data:
    input:
        # Path to the Python script responsible for scaling/normalizing the data
        script = "scripts/scaling_data.py",
        # Input file: the highly variable gene (HVG) filtered AnnData file
        hvg_h5ad = "data/merged_runs/merged_adata_hvg.h5ad",
    output:
        # Output file: the path for the normalized/scaled AnnData file
        scaled_h5ad = "data/merged_runs/merged_adata_scaled.h5ad"
    conda:
        # Specify the conda environment YAML file that contains the required dependencies
        "environment.yml"
    shell:
        # Command to execute the Python script, passing the input and output files as arguments
        "python {input.script} {input.hvg_h5ad} {output}"

# Rule to perform PCA on the normalized AnnData (H5AD) file
rule pca:
    input:
        # Path to the Python script responsible for performing PCA
        script="scripts/pca.py",
        # Input file: the normalized AnnData file
        scaled_h5ad = "data/merged_runs/merged_adata_scaled.h5ad"
    output:
        # Output path for the PCA-processed AnnData file
        pca_h5ad = "data/merged_runs/merged_adata_pca.h5ad"
        
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {input.scaled_h5ad} {output}"


# Rule to perform clustering and generate a UMAP plot
rule clustering_and_umap:
    input:
        # Path to the Python script responsible for clustering and UMAP
        script="scripts/clustering_and_umap.py",
        # Input file: the PCA-processed AnnData file
        pca_h5ad = "data/merged_runs/merged_adata_pca.h5ad",
        # Input: annotated AnnData file
        annotate_h5ad = "data/merged_runs/merged_PBMC_cell_type_annotations.csv"
    output:
        # Output paths for the UMAP plot and clustering results
        umap_h5ad = "data/merged_runs/merged_adata_umap.h5ad",
        umap_plot = "data/merged_runs/figures/umap_merged.png"
    params:
        umap_suffix = "_merged.png",
        path_to_save=config["working_directory"] + "/data/merged_runs/"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Execute the Python script with the input and output files
        "python {input.script} {input.pca_h5ad} {input.annotate_h5ad} {output.umap_h5ad} {output.umap_plot}"


# Rule to perform t-SNE clustering
rule t_sne:
    input:
        # Path to the Python script responsible for t-SNE clustering and generating the UMAP plot
        script = "scripts/t_sne.py",
        # Input file: the PCA-reduced AnnData file from previous preprocessing steps
        pca_h5ad = "data/merged_runs/merged_adata_pca.h5ad",
        # Input: annotated AnnData file
        annotate_h5ad = "data/merged_runs/merged_PBMC_cell_type_annotations.csv"
    output:
        # Output files: the AnnData file after t-SNE processing and the generated plot in PNG format
        t_sne_h5ad = "data/merged_runs/merged_adata_t_sne.h5ad",
        t_sne_plot = "data/merged_runs/figures/t_sne_merged.png"
    params:
        # Parameter for specifying suffix used in plot naming
        umap_suffix = "_merged.png",
        # Path to save output files, extracted from configuration
        path_to_save = config["working_directory"] + "/data/merged_runs/"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Run the Python script with the input PCA-reduced AnnData file and generate the t-SNE processed file and plot
        "python {input.script} {input.pca_h5ad} {input.annotate_h5ad} {output.t_sne_h5ad} {output.t_sne_plot}"


# Rule to perform differential expression analysis
rule differential_expression:
    input:
        # Path to the Python script that handles differential expression analysis, clustering, and UMAP visualization
        script="scripts/differential_expression.py",
        # Input file: the PCA-processed AnnData (h5ad) file, containing merged and dimensionally reduced data
        pca_h5ad = "data/merged_runs/merged_adata_pca.h5ad",
        # Input: annotated AnnData file
        annotate_h5ad = "data/merged_runs/merged_PBMC_cell_type_annotations.csv"
    output:
        # Output file path for the generated UMAP plot showing ranked gene groups after Leiden clustering
        de = "data/merged_runs/figures/rank_genes_groups_cell_type_merged.png"
    params:
        # Suffix used for the output differential expression plot
        de_suffix = "_merged.png",
        # Directory path where the results will be saved, taken from the configuration file
        path_to_save=config["working_directory"] + "/data/merged_runs/"
    conda:
        # Specify the conda environment for dependencies
        "environment.yml"
    shell:
        # Run the Python script, passing in the input PCA file, output suffix, and save directory as arguments
        "python {input.script} {input.pca_h5ad} {input.annotate_h5ad} {params.de_suffix} {params.path_to_save}"
