import scanpy as sc
import sys
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

warnings.filterwarnings("ignore")

#input_file = '/Users/emanuelmfonseca/project/single-cell-sequencing-migration/data/merged_runs/merged_adata_pca.h5ad'
#annotation_file = '/Users/emanuelmfonseca/project/single-cell-sequencing-migration/data/merged_runs/merged_PBMC_cell_type_annotations.csv'

def main(input_file, annotation_file, output_file, plot_output_file):
    # Load the preprocessed AnnData object
    adata = sc.read_h5ad(input_file)
    
    # Compute the neighborhood graph
    sc.pp.neighbors(adata)

    # Perform t-SNE analysis
    sc.tl.tsne(adata, n_pcs=50)  # Use the first 50 PCs for t-SNE

    # Apply clustering (Leiden)
    sc.tl.leiden(adata)

    # Load the cell type annotations from the saved CSV file
    annotations = pd.read_csv(annotation_file)

    # Add the annotations back to the AnnData object
    adata.obs['cell_type'] = annotations['cell_type'].to_list()

    # Save the filtered AnnData object with only highly variable genes
    adata.write(output_file)

   # Get the unique cell types
    cell_types = adata.obs['cell_type'].unique()
    
    # Generate a colormap for the cell types
    cmap = sns.color_palette("Set3", len(cell_types))
    
    # Create a mapping between cell types and the colors
    colors = {cell_type: cmap[i] for i, cell_type in enumerate(cell_types)} 
    
    # Create a t-SNE plot colored by clusters
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    for cell_type in cell_types:
        # Select the indices for the current cell type
        indices = adata.obs['cell_type'] == cell_type
        
        # Plot the points for this cell type
        ax.scatter(
            adata.obsm['X_tsne'][indices, 0],  # UMAP 1
            adata.obsm['X_tsne'][indices, 1],  # UMAP 2
            c=[colors[cell_type]], s=10, label=cell_type  # Use the mapped color
        )

    # Add labels and title for the plot
    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')
    ax.set_title('t-SNE Visualization of Clusters')
    
    # Add the discrete legend with cell type names and smaller font size
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cell Types', prop={'size': 5})

    # Adjust the layout
    plt.tight_layout()
    
    plt.savefig(plot_output_file, dpi=300)

if __name__ == "__main__":
    input_file = sys.argv[1]
    annotation_file = sys.argv[2]
    output_file = sys.argv[3]
    plot_output_file = sys.argv[4]
    main(input_file, annotation_file, output_file, plot_output_file)
