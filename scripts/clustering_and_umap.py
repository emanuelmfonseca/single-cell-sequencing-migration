import scanpy as sc
import sys
import warnings
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


warnings.filterwarnings("ignore")

input_file = "data/merged_runs/merged_adata_pca.h5ad"
annotation_file = "data/merged_runs/merged_PBMC_cell_type_annotations.csv"

def main(input_file, annotation_file, output_file, plot_output_file):
    # Read the input AnnData file
    adata = sc.read(input_file)
    
    # Compute the neighborhood graph of observations using PCA representation
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    # Compute UMAP (Uniform Manifold Approximation and Projection)
    sc.tl.umap(adata)
    
    # Perform Leiden clustering
    sc.tl.leiden(adata)
    
    # Load the cell type annotations from the saved CSV file
    annotations = pd.read_csv(annotation_file)

    # Add the annotations back to the AnnData object
    adata.obs['cell_type'] = annotations['cell_type'].to_list()

    # Write the clustering results to CSV files
    adata.write(output_file)

    # Get the unique cell types
    cell_types = adata.obs['cell_type'].unique()

    # Generate a colormap for the cell types
    cmap = sns.color_palette("Set3", len(cell_types))

    # Create a mapping between cell types and the colors
    colors = {cell_type: cmap[i] for i, cell_type in enumerate(cell_types)}

    # Create a UMAP plot
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot each cell type with its corresponding color
    for cell_type in cell_types:
        # Select the indices for the current cell type
        indices = adata.obs['cell_type'] == cell_type
        
        # Plot the points for this cell type
        ax.scatter(
            adata.obsm['X_umap'][indices, 0],  # UMAP 1
            adata.obsm['X_umap'][indices, 1],  # UMAP 2
            c=[colors[cell_type]], s=10, label=cell_type  # Use the mapped color
        )

    # Add labels and title for the plot
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('UMAP Visualization of Cell Types')

    # Add the discrete legend with cell type names and smaller font size
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cell Types', prop={'size': 5})

    # Adjust the layout to make room for the legend
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(plot_output_file, dpi=300)

if __name__ == "__main__":
    input_file = sys.argv[1]
    annotation_file = sys.argv[2]
    output_file = sys.argv[3]
    plot_output_file = sys.argv[4]
    main(input_file, annotation_file, output_file, plot_output_file)