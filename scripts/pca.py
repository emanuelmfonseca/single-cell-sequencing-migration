import scanpy as sc
import sys
import warnings

warnings.filterwarnings("ignore")

def main(input_file, output_file):
    # Read the input AnnData file
    adata = sc.read(input_file)
    
    # Perform PCA on the AnnData object
    sc.pp.pca(adata)
    
    # Save the processed AnnData object to the output file
    adata.write(output_file)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
