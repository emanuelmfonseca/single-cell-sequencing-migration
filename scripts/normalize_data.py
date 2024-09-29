import scanpy as sc
import sys
import warnings

warnings.filterwarnings("ignore")
        
def main(input_file, output_file):
    # Read the input AnnData file
    adata = sc.read(input_file)
    
    # Normalize total counts per cell to a target sum of 1e4
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Logarithmize the data
    sc.pp.log1p(adata)
    
    # Save the processed AnnData object to the output file
    adata.write(output_file)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
