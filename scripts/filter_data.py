import scanpy as sc
import sys
import warnings

warnings.filterwarnings("ignore")

def main(input_file,output_file):
    # Read the input file
    adata = sc.read(input_file)
    
    # Calculate metrics
    adata.obs['nCount_RNA'] = adata.X.sum(axis=1).A1  # Total counts per cell
    adata.obs['nFeature_RNA'] = (adata.X > 0).sum(axis=1).A1  # Number of genes per cell
    
    # Assuming mitochondrial genes are labeled with "MT-" in their names
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Identify mitochondrial genes
    adata.obs['percent.mt'] = adata[:, adata.var['mt']].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1 * 100
    
    # Define thresholds for filtering
    min_genes = 200      # Minimum number of detected genes (nFeature_RNA)
    max_genes = 2500     # Maximum number of detected genes (nFeature_RNA)
    max_counts = 25000   # Maximum number of total counts (nCount_RNA)
    max_mitochondria = 10  # Maximum percentage of mitochondrial genes (percent.mt)
    
    # Apply the filtering
    adata_filtered = adata[
        (adata.obs['nFeature_RNA'] > min_genes) &
        (adata.obs['nFeature_RNA'] < max_genes) &
        (adata.obs['nCount_RNA'] < max_counts) &
        (adata.obs['percent.mt'] < max_mitochondria)
    ].copy()
    
    # Save the preprocessed data to the output file
    adata_filtered.write(output_file)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)