import scanpy as sc
import sys
import warnings

warnings.filterwarnings("ignore")
        
def main(input_file, output_file):
    # Read the input AnnData file
    adata = sc.read(input_file)
    
    # Scale the data (mean 0, variance 1)
    sc.pp.scale(adata, max_value=10)
        
    # Save the processed AnnData object to the output file
    adata.write(output_file)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
