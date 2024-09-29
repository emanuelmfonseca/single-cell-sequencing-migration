import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main(input_file, output_file):
    # Load your AnnData object from Scanpy
    adata = sc.read_h5ad(input_file)

    # Calculate Total Counts and Number of Genes
    adata.obs['Total Counts'] = adata.X.sum(axis=1).A1  # Total counts per cell
    adata.obs['Number of Genes'] = (adata.X > 0).sum(axis=1).A1  # Number of genes per cell

    # Calculate Mitochondrial % (percentage of mitochondrial gene counts)
    # Assuming mitochondrial genes are labeled with "MT-" in their names
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Identify mitochondrial genes
    adata.obs['Mitochondrial %'] = adata[:, adata.var['mt']].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1 * 100

    # Create a figure for the merged plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Add letter labels for subplots
    letters = ['a)', 'b)', 'c)', 'd)', 'e)', '']
    for i, ax in enumerate(axes.flat):
        ax.text(-0.1, 1.1, letters[i], transform=ax.transAxes, size=20, weight='bold')

    # 1. Violin plots for Mitochondrial %, Number of Genes, and Total Counts
    sns.violinplot(y=adata.obs['Number of Genes'], ax=axes[0, 0], color="skyblue")
    axes[0, 0].set_title("Number of Genes")
    axes[0, 0].set_ylabel("Number of Genes")

    sns.violinplot(y=adata.obs['Total Counts'], ax=axes[0, 1], color="lightgreen")
    axes[0, 1].set_title("Total Counts")
    axes[0, 1].set_ylabel("Total Counts")

    sns.violinplot(y=adata.obs['Mitochondrial %'], ax=axes[0, 2], color="salmon")
    axes[0, 2].set_title("Mitochondrial %")

    # 2. Scatter plot: Number of Genes vs. Total Counts
    axes[1, 0].scatter(adata.obs['Total Counts'], adata.obs['Number of Genes'], c="blue", s=10, alpha=0.5)
    axes[1, 0].set_title("Number of Genes vs. Total Counts")
    axes[1, 0].set_xlabel("Total Counts")
    axes[1, 0].set_ylabel("Number of Genes")

    # 3. Scatter plot: Mitochondrial % vs. Total Counts
    axes[1, 1].scatter(adata.obs['Total Counts'], adata.obs['Mitochondrial %'], c="green", s=10, alpha=0.5)
    axes[1, 1].set_title("Mitochondrial % vs. Total Counts")
    axes[1, 1].set_xlabel("Total Counts")
    axes[1, 1].set_ylabel("Mitochondrial %")

    # Leave the last subplot (axes[1, 2]) empty
    axes[1, 2].axis('off')

    # Adjust the layout to prevent overlap
    plt.tight_layout()

    # Save the merged plot as an image
    plt.savefig(output_file, dpi=300)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)


