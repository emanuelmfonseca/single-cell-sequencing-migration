import scanpy as sc
import celltypist
import sys


def main(input_file, output_file):
    # Load your AnnData object from Scanpy
    adata = sc.read_h5ad(input_file)
    
    # Normalize total counts per cell to 10,000 reads
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Logarithmic transformation of the data
    sc.pp.log1p(adata)
    
    # Run CellTypist for PBMC cell type annotation using the Immune_All_Low model
    predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl')  # Model for PBMC annotation
    
    # Add the predicted cell type labels to the AnnData object
    adata.obs['cell_type'] = predictions.predicted_labels

    # Export the cell type annotations to a CSV file
    # This will save the cell type labels for each cell
    adata.obs[['cell_type']].to_csv(output_file)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)