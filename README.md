# Single-Cell Sequencing Migration Pipeline

This repository contains a bioinformatics pipeline designed to process single-cell sequencing data. The pipeline automates the processing of FASTQ files and manages the data using a cloud infrastructure. It is optimized for scalability and efficiency, reducing the need for manual intervention.

## Table of Contents
- [Introduction](#introduction)
- [Single-Cell Pipeline Overview](#single-cell-pipeline-overview)
- [AWS Architecture Overview](#aws-architecture-overview)
- [Cost Estimate](#cost-estimate)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
- [Main Findings](#main-findings)
- [Development and Testing Environment](#development-and-testing-environment)
- [Contributing](#contributing)
- [License](#license)

## Introduction
This pipeline automates the migration and processing of single-cell RNA sequencing data, providing a scalable and efficient solution for handling large datasets. Built with Snakemake and designed to run on AWS infrastructure, it leverages cloud services to streamline data processing, reduces manual intervention, and optimizes performance. 

The pipeline processes 10x Genomics data using Cell Ranger for primary analysis and Scanpy for secondary analysis, including quality control, normalization, dimensionality reduction, clustering, and differential expression. It supports seamless integration with AWS services like S3 for data storage and AWS Batch for orchestration, enabling fast, parallel processing of single-cell data.

This solution is ideal for automating large-scale single-cell sequencing projects while ensuring reliable and reproducible results.

## Single-Cell Pipeline Overview
This pipeline is designed to process single-cell RNA sequencing (scRNA-seq) data, performing both primary (blue) and secondary (orange) analyses. Key tasks include alignment, quality control, and advanced downstream analysis like clustering and differential expression.

![single-cell-pipeline](https://github.com/user-attachments/assets/f9e2c479-d0b6-457b-9fe6-8453e942e540)

### Primary Analysis

The primary analysis focuses on preparing raw single-cell RNA sequencing data for downstream analysis. This involves essential steps like data extraction, quality control, and basic processing, ensuring clean and structured datasets for secondary analysis.

- **Download Raw Data**
    - Raw sequencing data (FASTQ) is downloaded from the provided source for each run. These files are necessary for primary processing.

- **Reference Genome Download**
    - The human reference genome is downloaded and prepared for use in alignment and expression matrix generation during Cell Ranger processing.

- **Preprocessing with Cell Ranger**
    - The Cell Ranger tool is used to process the raw FASTQ files, aligning them to a reference genome, and generating gene expression matrices for each sample.

- **Quality Control Metrics**
    - Cells are evaluated based on several key indicators to ensure data quality. Cells with too few genes detected likely result from insufficient sequencing depth, while those with too many genes may represent doublets or multiplets. Additionally, cells with high mitochondrial RNA content are flagged as potentially stressed or degraded, optimizing the dataset for reliable downstream analysis.

- **Load and Process Data**
    - Process the Cell Ranger output to create AnnData (H5AD) objects. This step prepares the data for advanced analysis.

### Secondary Analysis

The secondary analysis builds on the results of the primary phase, involving deeper processing, normalization, scaling, and visualization of the single-cell RNA sequencing data. This stage ensures the data is ready for biological interpretation and further exploration.

- **Data Filtering**
    - Low-quality cells and genes are filtered out based on predefined criteria such as mitochondrial content and read counts, improving the dataset’s integrity.

- **Data Normalization**
    - Normalize the gene expression data to account for differences in sequencing depth across cells, preparing it for downstream analysis like clustering.

- **Feature Selection**
    - Identify highly variable genes across cells, focusing the analysis on genes that capture the most variance in the data. This helps to reduce noise and enhance the identification of key biological patterns.

- **Data Scaling**
    - Scale the gene expression data for each gene to have a mean of 0 and a standard deviation of 1. This step ensures that genes with different absolute expression levels contribute equally to downstream analyses, such as clustering and dimensionality reduction.

- **Data Annotation**
    - This step assigns cell type annotations to the filtered dataset using a custom Python script. It provides essential biological context, enabling downstream analyses such as clustering and differential expression.

- **Dimensionality Reduction (PCA)**
    - Perform Principal Component Analysis (PCA) on the filtered and normalized data to reduce dimensionality, facilitating visualization and clustering.

- **Clustering and Visualization**
    - Clustering is performed to group cells into distinct clusters based on gene expression patterns, and UMAP or t-SNE plots are generated for visualization.

- **Differential Expression Analysis**
    - Identify differentially expressed genes between cell clusters, providing insights into biological differences between cell populations.

## AWS Architecture Overview

### **Single-Cell Data Processing Pipeline Using AWS with IAM and Quilt Integration**:
This single-cell pipeline is designed to automate sequencing data processing using AWS services, Quilt for data management, and secure access management through IAM. The workflow covers data ingestion, processing, and storage, managed entirely through AWS infrastructure with Quilt providing version control and data lineage.

![aws-architecture](https://github.com/user-attachments/assets/6f6e4b11-f1b5-4628-abb2-92e1a7b7281f)

#### **1. Data Ingestion into S3 and Quilt**:
- Sequencing data is uploaded to an **S3 bucket** and registered in **Quilt** for data management. S3 serves as the central storage for raw sequencing data (e.g., FASTQ files), intermediate results, and final outputs, while Quilt ensures metadata tracking and versioning.
- **IAM roles** control access to S3 and Quilt, allowing only authorized services and users to read and write data.
- Every time new data is uploaded, an event triggers AWS services to start the pipeline.

#### **2. Lambda Trigger**:
- **AWS Lambda** is automatically triggered when new data is uploaded to **S3** and registered in Quilt. Lambda functions, using **IAM roles**, have the permissions required to interact with both **S3** and Quilt, execute workflows, and trigger downstream processes (e.g., EC2 or AWS Batch).
- This automation ensures the pipeline starts processing new data as soon as it becomes available.

#### **3. Running the Snakemake Pipeline on EC2**:
- **Lambda**, using its designated **IAM role**, triggers the launch of an **EC2 instance** to execute the Snakemake pipeline. The **EC2 instance** is assigned an **IAM role** with permissions to access **S3** and **Quilt** for retrieving raw data, running the pipeline, and sending logs and metrics to **CloudWatch**.
- The pipeline, running on EC2, downloads raw data from S3 via Quilt, processes it (e.g., quality control and alignment), and tracks metadata updates.

#### **4. Storing Outputs in S3 and Quilt**:
- After processing, output files (e.g., BAM, H5AD files) are uploaded back to **S3** and registered in **Quilt** for versioning and data management.
- The **EC2 instance**, using its **IAM role**, ensures that results are securely stored in both **S3** and Quilt for traceability.

#### **5. Batch Processing**:
- For larger datasets, **AWS Step Functions** will be used to orchestrate parallel processing. Step Functions will manage the workflow by submitting jobs to **AWS Batch**, which provisions **EC2** instances to run the Snakemake pipeline. **IAM roles** assigned to the **AWS Batch** instances will ensure secure access to **S3** data, interaction with **Quilt** for metadata management, and the execution of the Snakemake workflows.

#### **6. Workflow Orchestration with Step Functions**:
- **AWS Step Functions** orchestrate the pipeline, ensuring tasks such as data retrieval, processing, and storage are completed sequentially.
- **IAM roles** allow Step Functions to manage permissions between Lambda, EC2, Quilt, and other services, enabling secure and orderly workflow execution.

#### **7. Monitoring and Logging with CloudWatch**:
- Logs from EC2, Lambda, and Batch instances are sent to **CloudWatch** for real-time monitoring.
- **IAM roles** ensure secure transmission of logs and metrics, enabling detailed monitoring and troubleshooting.

#### **8. Secure Access with IAM**:
- **IAM roles** govern access between AWS services and Quilt, ensuring that each service (Lambda, EC2, Batch) has the necessary permissions to interact with other AWS resources securely.
- **IAM roles** are assigned across the pipeline, providing granular control over access to S3, Quilt, Lambda, EC2, Step Functions, and CloudWatch to ensure robust security. 

### Updated **Data Flow Summary with Quilt**:
1. **S3 → Quilt → Lambda**: Data is uploaded to S3 and registered in Quilt. This triggers a Lambda function, with **IAM roles** granting access to both S3 and Quilt.
2. **Lambda → EC2**: Lambda starts an EC2 instance, and **IAM roles** allow the EC2 to fetch data from Quilt (via S3) and run the Snakemake pipeline.
3. **EC2 → Quilt → S3**: Processed results are registered back into Quilt and stored in S3 using the **IAM role** assigned to EC2.
4. **EC2 → CloudWatch**: Logs from EC2 are sent to CloudWatch using **IAM permissions**.
5. **AWS Batch** uses **IAM roles** to manage job scheduling for parallel processing.
6. **Step Functions → Lambda/EC2**: Step Functions orchestrate the workflow using **IAM roles** for managing permissions between all services.

## Cost Estimate

**Here are preliminary cost estimates for the key AWS services involved in the single-cell pipeline, as calculated using the AWS Pricing Calculator:**

| Service                    | Details                                                                                                         | Price per Run ($) | Value per Month ($) |
|----------------------------|-----------------------------------------------------------------------------------------------------------------|-------------------|---------------------|
| **Amazon S3 (Storage)**     | - **Storage**: 1.5 TB of data per month <br> - **Price for S3 Standard**: ~$0.023 per GB/month <br> - **Intermediate files**: FASTQ, BAM, H5AD files are stored in intermediate stages | 3.24              | 103.68              |
| **Amazon S3 (Retrieval)**   | - **Data Retrieval**: 1.5 TB retrieved per month <br> - **Standard retrieval price**: ~$0.01 per GB <br> - **Retrieval cost**: $15.36 (1.5 TB * 1024 GB/TB * $0.01) | 0.48              | 15.36               |
| **EC2 (Snakemake pipeline)**| - **Instance type**: m5.xlarge (4 vCPUs, 16 GB RAM) <br> - **Estimated hours per sequencing run**: 36 hours per run <br> - **Runs per week**: 8 runs (total 32 runs/month) <br> - **Total hours per month**: 36 hours per run * 32 runs/month = 1,152 hours/month <br> - **Price for m5.xlarge**: ~$0.192 per hour <br> - **Monthly EC2 cost**: 1,152 hours * $0.192 = $221.18 | 6.91              | 221.18              |
| **Lambda (Trigger function)**| - **Number of invocations**: 32 per month (1 per sequencing run) <br> - **Execution time**: ~1 second (1000 ms) <br> - **Memory**: 128 MB <br> - **Price**: $0.00001667 per request (100 ms increments) <br> - **Monthly Lambda cost**: $0.05 (32 invocations * $0.00001667) | 0.0019            | 0.05                |
| **AWS Batch**              | - **Assume 10 EC2 instances (m5.xlarge)** running in parallel for batch jobs <br> - **Total runtime for batch jobs**: 12 hours per run <br> - **Total EC2 hours for Batch**: 10 instances * 12 hours = 120 hours per run <br> - **For 32 runs per month**: 120 hours per run * 32 runs = 3,840 hours per month <br> - **Price for m5.xlarge**: ~$0.192 per hour <br> - **Monthly Batch cost**: 3,840 hours * $0.192 = $737.28 | 23.04             | 737.28              |
| **CloudWatch**             | - **Logs volume**: 15 GB/month <br> - **Custom metrics**: 10 metrics/month <br> - **Log storage cost**: ~$0.50 per GB <br> - **Metrics cost**: ~$0.30 per metric <br> - **Monthly CloudWatch cost**: $12.50 (15 GB * $0.50 + 10 metrics * $0.30) | 0.78              | 12.50               |
| **Step Functions**         | - **State transitions**: 32 transitions (1 per sequencing run) <br> - **Cost per transition**: $0.025 per 1,000 transitions <br> - **Monthly Step Functions cost**: $0.02 (32 transitions * $0.025/1000) | 0.0004            | 0.02                |

**Total Cost per Run**: $33.71  
**Total Monthly Cost**: $1,089.07

## Installation

### Dependencies

To run this project, you will need the following dependencies installed:

- **Conda**: An environment manager to handle package installation and dependency management.

**Note**: All other dependencies are specified in the `yaml` file and will be automatically installed when the environment is created, ensuring that all necessary tools and libraries are properly set up.

### Installation Instructions

**Follow these steps to install Conda and access the pipeline:**

1. Ensure you have Conda installed. If not, you can download it from [here](https://docs.conda.io/en/latest/miniconda.html).

2. Clone the repository and navigate to the project directory:

    ```bash
    git clone https://github.com/emanuelmfonseca/single-cell-sequencing-migration.git
    cd single-cell-sequencing-migration
    ```

3. Create and activate the Conda environment using the provided `environment.yml` file:

    ```bash
    conda env create -f environment.yml
    conda activate single-cell-sequencing-migration
    ```

This will automatically install **Snakemake** and any other necessary packages specified in the environment configuration. **Snakemake** is a powerful workflow management system designed to create reproducible, scalable, and automated data analysis pipelines.


4. **Install Cell Ranger**:  
   
   To process the single-cell RNA sequencing dataset, you will need to install **Cell Ranger**, a key tool developed by 10x Genomics for analyzing single-cell data. Please follow the detailed installation instructions provided by 10x Genomics:

   **[Cell Ranger Installation Guide](https://www.10xgenomics.com/support/software/cell-ranger/latest)**  

   Ensure that **Cell Ranger** is correctly installed and added to your system's PATH so that it can be accessed from the command line within your pipeline. For more information on configuring PATH variables, refer to the relevant documentation on your operating system.

## Usage Guide

### 1. Data Preparation

The input dataset for this pipeline will be automatically downloaded from 10x Genomics website. Specifically, the pipeline retrieves the **[1k PBMCs from a Healthy Donor (v3 Chemistry)](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)** dataset. No manual download is required, as the pipeline includes a step to fetch this dataset during execution.


### 2. Configuration

Before running the pipeline, you must configure the `config.yaml` file. This file contains essential parameters that control how the pipeline executes.

#### Key Configuration Parameters:

- **`working_directory`**: Full path to your project’s working directory.
- **`weekly_runs`**: Number of independent genomic runs to be performed.
- **`cell_ranger_path`**: Full path to the Cell Ranger installation directory. This is required for processing single-cell RNA sequencing data, ensuring compatibility with the 10x Genomics pipeline.
- **`ncores`**: Number of CPU cores to be used for parallel processing.

**Note**: Make sure to correctly set the `working_directory` to avoid pipeline failures.

### 3. Running the Pipeline

Once the configuration is complete, you can execute the pipeline using Snakemake.

#### Command to Run:

```bash
snakemake --cores <number_of_cores>
```

#### Example:

If you want to run the pipeline using 6 CPU cores:

```bash
snakemake --cores 6
```

**Tip**: Adjust the number of cores based on your system’s capacity to optimize performance.

### 4. Troubleshooting & Tips

- Ensure all required software, including Cell Ranger, is installed before running the pipeline.
- Double-check the file paths and parameters in `config.yaml` to prevent issues related to misconfiguration.
- Test the pipeline with a smaller dataset or fewer cores to validate the setup before full execution. Use Snakemake’s --dry-run to check dependencies and rule definitions without running the actual commands.

### Outputs

- The pipeline generates various outputs at each step, including:

  - **Cell Ranger output files**: Processed single-cell sequencing data using Cell Ranger, including raw count matrices.
    - Example: `data/run_1/cellranger/outs`.

  - **Reference genome files**: Human reference genome used for alignment and analysis.
    - Example: `data/human_reference/`.

  - **Quality Control (QC) plots**: Visualizations of unfiltered and filtered QC metrics such as mitochondrial percentage and feature counts.
    - Example: `data/merged_runs/figures/qc_metrics_plot_unfiltered.png`.
    - Example: `data/merged_runs/figures/qc_metrics_plot_filtered.png`.

  - **Normalized AnnData (H5AD) files**: Processed single-cell RNA-seq data after normalization.
    - Example: `data/merged_runs/merged_adata_normalized.h5ad`.

  - **Annotated AnnData (H5AD) files**: Annotated single-cell RNA-seq data based on cell types.
    - Example: `data/merged_runs/merged_adata_annotate_h5ad.h5ad`.
    - Example: `data/merged_runs/merged_PBMC_cell_type_annotations.csv`.

  - **Clustered AnnData files**: Files containing clustering information for the cells.
    - Example: `data/merged_runs/merged_adata_pca.h5ad`.

  - **UMAP plot**: UMAP visualization of cell clusters based on transcriptomic data.
    - Example: `data/merged_runs/figures/umap_merged.png`.

  - **t-SNE plot**: t-SNE visualization of cell clusters for comparison with UMAP results.
    - Example: `data/merged_runs/figures/t_sne_merged.png`.

  - **Differential Expression (DE) plots**: Visualizations showing ranked gene groups after clustering.
    - Example: `data/merged_runs/figures/rank_genes_groups_cell_type_merged.png`.


## Main Findings

### Unfiltered QC Metrics
1. **Number of Genes**: A broader distribution in unfiltered data shows more low-quality cells with fewer genes detected.
2. **Total Counts**: Higher variability in total counts in unfiltered cells suggests the presence of poor-quality cells.
3. **Mitochondrial %**: Unfiltered data includes cells with high mitochondrial expression, indicating potential cell stress or damage.

- Unfiltered Data shows violin and scatter plots summarizing quality control metrics before filtering. Panel **a** presents gene count distribution, while panel **b** shows total RNA counts per cell. Panel **c** illustrates mitochondrial gene percentages. Panel **d** shows the correlation between gene counts and total RNA counts, and panel **e** presents the relationship between mitochondrial percentage and total RNA counts.

![qc_metrics_plot_unfiltered](https://github.com/user-attachments/assets/b40d3a51-ca07-487a-aaf0-d0a606d47715)

### Filtered QC Metrics
1. **Number of Genes**: The filtered data shows a tighter range of gene counts per cell, centering around a higher median, indicating higher-quality cells.
2. **Total Counts**: Filtered cells display higher RNA counts and a more consistent distribution, removing outliers from low-quality cells.
3. **Mitochondrial %**: Lower mitochondrial percentages in filtered cells suggest removal of stressed or dying cells.

- Filtered Data displays similar plots after filtering. Panel **a** presents the gene count distribution, panel **b** shows total RNA counts, and panel **c** illustrates mitochondrial percentages. Panel **d** shows the correlation between gene counts and total counts, while panel **e** presents the relationship between mitochondrial percentage and total RNA counts.

![qc_metrics_plot_filtered](https://github.com/user-attachments/assets/83895528-1cee-4641-8b89-f5b138127e98)
 
### UMAP Clustering
- UMAP emphasizes the global structure of the data, showing transitions between cell types. Cells are color-coded by type, and clusters indicate distinct but related populations. UMAP can highlight gradual changes, such as cell differentiation pathways.

 - UMAP presents a UMAP projection, highlighting cellular diversity and clustering by cell type.

![umap_merged](https://github.com/user-attachments/assets/14c46ea5-1b4c-473e-b7dc-370c2f65ea77)

### t-SNE Clustering
- Cells cluster based on gene expression profiles, with distinct groups representing different immune cell types. t-SNE shows clear separation between clusters, suggesting high-resolution identification of individual populations in the dataset.

- t-SNE shows a t-SNE plot clustering cells by type, providing a complementary view to the UMAP.

![t_sne_merged](https://github.com/user-attachments/assets/3e871ebd-b6af-4d61-97f9-c65a70e9f60d)

### Rank Genes by Cell Type
- **Top Genes**: Each plot shows genes ranked by significance for identifying specific cell types. Top marker genes for each population indicate distinct expression profiles, providing clear differentiation between cell types (e.g., CD16+ NK Cells vs. other cells).

- Ranked Genes displays ranked gene expression differences for each cell type compared to all other cells.

![rank_genes_groups_cell_type_merged](https://github.com/user-attachments/assets/98360dc5-8fbd-4404-9566-96dc0abc2c2f)

## Development and Testing Environment
This pipeline was developed on a MacBook 2020 with an M1 chip. However, the Cell Ranger part of the pipeline was tested and run exclusively on a Linux machine, as Cell Ranger is only supported on Linux systems. The tutorial associated with this project was tested on a MacBook with the same configuration as the development environment.

## Contributing
Contributions are welcome. Please open an issue or submit a pull request for any improvements. For questions, feel free to email emanuelmfonseca@gmail.com.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.
