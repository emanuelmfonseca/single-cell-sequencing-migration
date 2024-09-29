import scanpy as sc
import sys
import warnings

warnings.filterwarnings("ignore")

def main(output_file, *input_files):
    # Read all the input AnnData files
    adata_files = [sc.read(file) for file in input_files]
    
    # Create batch keys based on the number of input files
    runs = [f"run{run_id}" for run_id in range(1, len(adata_files) + 1)]

    # Merge them with Scanpy's concat function
    merged_data = sc.concat(adata_files, label='batch', keys=runs)

    # Save the merged AnnData object
    merged_data.write(output_file)

if __name__ == "__main__":
    output_file = sys.argv[1]
    input_files = sys.argv[2:]  # Collect all input files
    main(output_file, *input_files)
