import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from visualisation_functions import *

def calculate_pca_browning_score(protein_df, n_components=2):
    scores = pd.DataFrame(index=protein_df.index)
    print(protein_df.shape)
    protein_df.drop(columns=['Batch'], inplace=True, errors='ignore')  # Ensure 'Batch' column is dropped if present
    scaler = StandardScaler()
    protein_df_scaled = scaler.fit_transform(protein_df)
    
    pca = PCA(n_components=n_components, random_state=42) # Added random_state for reproducibility
    principal_components = pca.fit_transform(protein_df)

    pc_df = pd.DataFrame(data=principal_components, 
                         columns=[f'PC{i+1}' for i in range(n_components)],
                         index=protein_df.index)
    

    for i in range(n_components):
        scores[f'pca_browning_score_PC{i+1}'] = pc_df.iloc[:, i]
    #scores[f'pca_browning_score_PC{component_to_use+1}'] = pc_df.iloc[:, component_to_use]
    
    print(f"PCA Explained Variance Ratios: {pca.explained_variance_ratio_}")
    pca_feature_loadings = pd.DataFrame(
        pca.components_.T,
        index=protein_df.columns,
        columns=[f'PC{i+1}' for i in range(n_components)]
    )
    return scores, pca_feature_loadings#, pc_df, pca_feature_loadings


def calculate_composite_score(all_scores_df, weights=None):
    composite_score = pd.Series(0.0, index=all_scores_df.index)
    total_weight = 0.0 # Ensure float for potential division
    
    # Create a new series for calculation, ensure it's float
    # Initialize with 0.0 to handle cases where some samples might have all NaN scores for weighted items
    calculated_composite_score = pd.Series(0.0, index=all_scores_df.index, dtype=float)
    
    # Keep track of how many weighted scores contributed to each sample's composite score
    # This is useful if we want to normalize by sum of weights of *contributing* scores
    # For now, simpler: normalize by sum of all *intended* weights (total_weight)
    if weights == None:
        weights = {i:1 for i in all_scores_df.columns}

    for score_name, weight in weights.items():
        if score_name in all_scores_df.columns:
            # Fill NaN with 0 for summation; consider other imputation if 0 is not appropriate
            calculated_composite_score += all_scores_df[score_name].fillna(0) * weight 
            total_weight += weight
        else:
            print(f"Warning: Score '{score_name}' not found in DataFrame for composite score.")
            
    if total_weight == 0:
        print("Warning: Total weight is zero, composite score cannot be calculated meaningfully.")
        return pd.Series(np.nan, index=all_scores_df.index) # Return NaNs
    
    # Optional: Normalize by total_weight if weights don't sum to 1, or to average contributions
    # If a sample had all NaNs for the input scores, fillna(0) makes them 0.
    # If sum of weights is not 1, this normalization makes sense.
    # If sum of weights IS 1, this division isn't strictly necessary unless to handle precision.
    # For this version, let's assume weights are proportions summing to 1, or are relative.
    # If they don't sum to 1, the scale of composite_score depends on sum of weights.
    # For now, return the weighted sum directly. User can normalize later if needed.
    
    return pd.DataFrame(calculated_composite_score, index = calculated_composite_score.index, columns=['composite_score'])
