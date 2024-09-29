import scanpy as sc
import sys
import warnings

warnings.filterwarnings("ignore")

def main(input_file, output_file):
    # Load the 10X matrix using gene symbols for variable names and enable caching for faster access
    adata = sc.read_10x_mtx(input_file, var_names="gene_symbols", cache=True)

    # Write the loaded data into an H5AD file for further analysis
    adata.write(output_file)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
