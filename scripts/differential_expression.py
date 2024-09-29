import scanpy as sc
import sys
import warnings
import os
import pandas as pd

warnings.filterwarnings("ignore")

def main(input_file, annotation_file, output_file_de, working_directory):
    # Load the preprocessed AnnData object
    adata = sc.read_h5ad(input_file)
    
    # Perform logarithmic transformation
    sc.pp.log1p(adata)
    
    # Compute nearest neighbors
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Perform Leiden clustering
    sc.tl.leiden(adata)

    # Load the cell type annotations from the saved CSV file
    annotations = pd.read_csv(annotation_file)

    # Add the annotations back to the AnnData object
    adata.obs['cell_type'] = annotations['cell_type'].to_list()
    
    # Rank genes using Wilcoxon test for each cluster
    sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon')
    
    os.chdir(working_directory)

    # Save the plot of the top 10 ranked genes for each cluster
    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, save=output_file_de, show=False)

if __name__ == "__main__":
    input_file = sys.argv[1]
    annotation_file = sys.argv[2]
    output_file_de = sys.argv[3]
    working_directory = sys.argv[4]
    main(input_file, annotation_file, output_file_de, working_directory)

