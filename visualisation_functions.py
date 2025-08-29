import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE as Tsne
from sklearn.preprocessing import StandardScaler, LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression



def plot_pca_pcs(pca_scores, markers, sample_labels, batch_labels):
    """Plots PCA of the data with different markers and colors based on tissue and diet."""
    fig, axs = plt.subplots(ncols=2, figsize=(11,5))
    scale = StandardScaler()

    x1=sns.scatterplot(x=pca_scores.pca_browning_score_PC1, y=pca_scores.pca_browning_score_PC2, hue=sample_labels, ax = axs[0], legend=True)
    x2=sns.scatterplot(x=pca_scores.pca_browning_score_PC1, y=pca_scores.pca_browning_score_PC2, hue=batch_labels, ax = axs[1], legend=True)
    fig.suptitle("PCA")
    fig.tight_layout(rect=[0, 0, 0.9, 1]) # Adjust the right boundary of the tight_layout rectangle
    plt.savefig("script_outputs/pca_plot.png", dpi=300)
    print('Plot saved as script_outputs/pca_plot.png')
    plt.close()

def plot_score_distributions(df, status_labels):
    """Plots histograms and boxplots for each score, grouped by status if available."""
    print("\n## Plotting Score Distributions ##")
    for col in df.columns:
        plt.figure(figsize=(12, 5))
        
        # Histogram
        plt.subplot(1, 2, 1)
        sns.histplot(df[col].dropna(), kde=True)
        plt.title(f"Histogram of {col}")
        plt.xlabel(col)
        plt.ylabel("Frequency")
        
        # Boxplot by Status (if status_labels are available)
        if status_labels is not None:
            plt.subplot(1, 2, 2)
            # Combine scores and status for plotting
            plot_df = pd.DataFrame({'score': df[col], 'Status': status_labels})
            plot_df['score'] = pd.to_numeric(plot_df['score'], errors='coerce')
            plot_df = plot_df.dropna(subset=['score'])

            # Define order for consistent plotting
            status_order = ['White', 'Intermediate', 'Brown'] 
            # Filter out statuses not present to avoid errors with order
            present_statuses = [s for s in status_order if s in plot_df['Status'].unique()]

            if present_statuses:
                 sns.boxplot(x='Status', y='score', data=plot_df, order=present_statuses, palette="Set2", hue="Status")
                 #sns.stripplot(x='Status', y='score', data=plot_df, order=present_statuses, color=".25", size=3, alpha=0.5)
                 plt.title(f"Boxplot of {col} by Status")
                 plt.xlabel("Sample Status")
                 plt.ylabel(col)

                 # Perform Mann-Whitney U tests for Brown vs White
                 if 'Brown' in present_statuses and 'White' in present_statuses:
                     brown_scores = plot_df[plot_df['Status'] == 'Brown']['score'].dropna().values
                     white_scores = plot_df[plot_df['Status'] == 'White']['score'].dropna().values
                     if len(brown_scores) > 1 and len(white_scores) > 1:
                         try:
                             stat, p_val = mannwhitneyu(brown_scores, white_scores, alternative='two-sided')
                             plt.text(0.95, 0.05, f'Brown vs White P-val: {p_val:.2e}',
                                      ha='right', va='bottom', transform=plt.gca().transAxes, fontsize=9)
                         except ValueError as e:
                             print(f"Could not compute Mann-Whitney U for {col} (Brown vs White): {e}")


            else: # If no relevant statuses for boxplot
                plt.text(0.5, 0.5, "Status information not available\nor no relevant statuses for boxplot.",
                         ha='center', va='center', transform=plt.gca().transAxes)

        plt.tight_layout()
        plt.savefig('script_outputs/' + col + '_distribution.png', dpi=300)
        plt.close()
        

def get_base_plots(protein_df, sample_labels, batch_labels, tissue, diet, all=None, markers=None, lasso=None, scaling=True):
    """Generates base plots for the protein data, including PCA and UMAP visualizations."""
    print("\n## Generating Base Plots ##")
    dim_reducers = [UMAP(n_components=2, random_state=1), PCA(n_components=2), Tsne(n_components=2, random_state=1)]
    reducer_names = ['UMAP', 'PCA', 't-SNE']
    scaler = StandardScaler()
    
    for reducer, name in zip(dim_reducers, reducer_names):
        fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(16,16))
        for row, label in enumerate([sample_labels, batch_labels, diet]):
            for col, (gset, setname) in enumerate(zip([all, markers, lasso], ['All Genes', 'Markers', 'LASSO'])):
                ax = axs[row, col]
                #hue = sample_labels if isinstance(label, pd.Series) else label.values
                if scaling:
                    scaled_data = scaler.fit_transform(protein_df[gset]) if gset else scaler.fit_transform(protein_df)
                else:
                    scaled_data = protein_df[gset] if gset else protein_df
                reduced_data = reducer.fit_transform(scaled_data)
                sns.scatterplot(x=reduced_data[:, 0], y=reduced_data[:, 1], hue=label, ax=axs[row, col], legend=False if col < 2 else True)
                axs[row, col].set_title(f"{name} - {label.name} - {setname if gset else 'All Genes'}")
        # Add legend only for the rightmost plots
                if col == 2:
                    if row == 0:
                        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
                    elif row == 1:
                    # For detailed groups, show legend outside
                        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=6)
                    else:
                        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
        plt.tight_layout()
        plt.savefig(f'script_outputs/base_plots/{name}_base_plots_scaling{str(scaling)}.png', dpi=300)
        plt.close()

    
def test_umap(data, markers):
    metadata = pd.read_csv('../data/ml_training/meta107_samples.csv', index_col=0, header=0)
    ld = {s:l for s, l in zip(metadata.file_name, metadata["tissue"])}
    ldd = {s:l for s, l in zip(metadata.file_name, metadata["label"])}
    
    fig, axs = plt.subplots(ncols=3, figsize=(16,5))

    hue=[ldd[s] for s in data.index]
    pca = UMAP(n_components=2)
    pcad = pca.fit_transform(data)


    x=sns.scatterplot(x=pcad[:,0], y=pcad[:,1], hue=hue, ax = axs[0], legend=True)
    axs[0].set_title("All genes")
    hue=[ldd[s] for s in data.index]
    sns.scatterplot(x=pcad[:,0], y=pcad[:,1], hue=[ld[s] for s in data.index], ax = axs[1], legend=True,)
    hue=[ldd[s] for s in data.index]
    pcad = pca.fit_transform(data[markers])
    sns.scatterplot(x=pcad[:,0], y=pcad[:,1], hue=[ld[s] for s in data.index], ax = axs[2], legend=False, )
    axs[2].set_title("Marker genes")
    fig.suptitle("UMAP")

    # Get unique hue categories and assign colors
    unique_hues = sorted(list(set(hue)))
    palette = sns.color_palette(n_colors=len(unique_hues))
    hue_map = dict(zip(unique_hues, palette))

    # Create legend handles manually
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=label,
                        markerfacecolor=hue_map[label], markersize=8)
            for label in hue]

    # Add the legend to the figure, placing it to the right
    #fig.legend(handles=handles, title="Tissue-Diet", loc='center left', bbox_to_anchor=(0.9, 0.5))
    #fig.legend(handles=x.handles, labels=hue)
    # Adjust layout to prevent the legend overlapping the plots
    fig.tight_layout(rect=[0, 0, 0.9, 1]) # Adjust the right boundary of the tight_layout rectangle
    plt.show()