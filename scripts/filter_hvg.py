import scanpy as sc
import sys
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

def main(input_file, output_file, plot_output_file):
    # Load the preprocessed AnnData object
    adata = sc.read_h5ad(input_file)
    
    # Logarithmic transformation of the data
    sc.pp.log1p(adata)
    
    # Detect highly variable genes (HVGs)
    sc.pp.highly_variable_genes(adata, min_mean=0.001, max_mean=3, min_disp=0.5)
    
    # Save the filtered AnnData object with only highly variable genes
    adata[:, adata.var.highly_variable].write(output_file)

    # Create a figure for hvg
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))

    # Plot non-highly variable genes (set them to light gray for distinction)
    axes.scatter(
        adata.var['means'][~adata.var['highly_variable']], 
        adata.var['dispersions_norm'][~adata.var['highly_variable']], 
        color='lightgray', s=10, label='Non-highly variable genes'
    )

    # Plot highly variable genes
    axes.scatter(
        adata.var['means'][adata.var['highly_variable']], 
        adata.var['dispersions_norm'][adata.var['highly_variable']], 
        color='red', s=10, label='Highly variable genes'
    )

    # Add labels and title for the first panel
    axes.set_xlabel('Mean expression')
    axes.set_ylabel('Normalized dispersion')
    axes.set_title('Highly Variable Genes')

    # Add the legend
    axes.legend()

    # Adjust the layout
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(plot_output_file, dpi=300)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_output_file = sys.argv[3]
    main(input_file, output_file, plot_output_file)
