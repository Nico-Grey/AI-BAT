import pandas as pd
import numpy as np
def load_example_data(protein_data_path = '../data/data_python/corrected.csv', precomputed_lasso_coefs_path='../data/data_python/coefs_from_lasso.pkl', 
                      brown_markers_path='../data/data_python/marker.txt', sample_labels_path='../data/data_python/meta107_samples.csv'):
    
    protein_df = pd.read_csv(protein_data_path, index_col=0, header=0).T

    coefs = pd.read_pickle(precomputed_lasso_coefs_path).T

    brown_markers_from_file = pd.read_csv(brown_markers_path, sep="\t", index_col=0, header=0)["GeneName"].to_list()
    brown_markers_from_file = [i for i in brown_markers_from_file if i in protein_df.columns]
    

    metadata = pd.read_csv(sample_labels_path, index_col=0, header=0)
    batchd = {s:s.split('-')[0] for s in metadata.index}
    tissue = pd.Series(metadata.tissue, index=metadata.index, name="Tissue")
    diet = pd.Series(metadata.diet, index=metadata.index, name="Diet")
    batch = [batchd[s] for s in protein_df.index]

    wc = [i for j, i in zip(coefs.eWAT, coefs.index) if abs(j) > 0.06]
    bc = [i for j, i in zip(coefs.BAT, coefs.index) if abs(j) > 0.06]
    
    brown_markers_list = brown_markers_from_file + bc
    white_markers_list = wc
    
    sample_labels = pd.read_csv(sample_labels_path, index_col=1, header=0)
    d = {"BAT": "Brown", "eWAT": "White", "WAT": "White", "iWAT": "Intermediate", "sWAT":"Intermediate"}
    sl = [d[i] for i in sample_labels.tissue.to_list()]
    sample_labels_series = pd.Series(sl, index=sample_labels.index, name="Status")
    batch_labels_series = pd.Series(batch, index=sample_labels.index, name="Batch")
    lda_features = brown_markers_list + white_markers_list

    return protein_df, brown_markers_list, white_markers_list, sample_labels_series, batch_labels_series, lda_features, brown_markers_from_file, tissue, diet
