"""
AI-BATS Browning Pipeline - Machine Learning Analysis
====================================================

This script implements a comprehensive machine learning pipeline for analyzing
protein expression data related to browning and whitening of adipose tissue.

The pipeline performs:
1. Data loading and preprocessing
2. Feature selection using LASSO regression
3. Dimensionality reduction (PCA) for feature engineering
4. Multiple machine learning models (LDA, Random Forest, TabPFN, Logistic Regression)
5. Cross-validation and scoring
6. Composite score calculation
7. Visualization and analysis output

Key Models Used:
- Linear Discriminant Analysis (LDA): For dimensionality reduction and classification
- Random Forest: For ensemble-based classification with feature importance
- TabPFN: For state-of-the-art tabular data classification
- Logistic Regression: For interpretable linear classification
- PCA: For unsupervised dimensionality reduction

Output Files:
- all_sample_scores.csv: Comprehensive scoring results for all samples
- Various PNG plots: Feature importance, distributions, correlations
- SHAP analysis plots: Model interpretability visualizations
- Base plots: PCA, t-SNE, UMAP visualizations

Dependencies:
- Custom modules: data_functions, ml_functions, scoring_functions, 
                 visualisation_functions, analysis_functions
- External libraries: sklearn, pandas, numpy, matplotlib, seaborn, 
                     umap, shap, tabpfn, scipy

Author: AI-BATS Team
Version: 1.0
"""

# ==============================================================================
# IMPORT STATEMENTS
# ==============================================================================

# Core machine learning and data processing libraries
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis 
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, cross_val_score

# Data manipulation and numerical computing
import numpy as np
import pandas as pd

# Visualization libraries
import matplotlib.pyplot as plt
import seaborn as sns

# Advanced ML and analysis libraries
import umap  # Uniform Manifold Approximation and Projection
import shap  # SHapley Additive exPlanations for model interpretability

# Specialized ML models
from tabpfn import TabPFNClassifier  # Prior-Fitted Networks for tabular data
from scipy.stats import mannwhitneyu  # Statistical testing

# File handling
from pathlib import Path
import pickle

# Custom AI-BATS modules containing specialized functions
from data_functions import *           # Data loading and preprocessing utilities
from ml_functions import *            # Machine learning model implementations
from scoring_functions import *       # Scoring and evaluation functions
from visualisation_functions import * # Plotting and visualization utilities
from analysis_functions import *      # Statistical analysis functions

# ==============================================================================
# CONFIGURATION AND PARAMETERS
# ==============================================================================

# This script is designed to run a complete analysis pipeline for protein expression data related to browning and whitening of adipose tissue.

# Example weights for composite score calculation
# These weights determine how different scoring methods contribute to the final composite score
# Positive weights favor browning, negative weights favor whitening
example_weights = {
        'PC2_markers': 0.1,              # PCA Component 2 based on marker genes
        'PC2': -0.1,                     # PCA Component 2 based on all features
        #'PC1_markers': 0.1,             # PCA Component 1 based on marker genes (commented out)
        #'PC1': 0.1,                     # PCA Component 1 based on all features (commented out)
        'LogisticRegression_probability_score': 0.4,  # Logistic regression probability
        'RandomForest_probability_score': 0.5,        # Random forest probability
        'tabpfn_score_Brown': 3,         # TabPFN score for brown classification (high weight)
        'tabpfn_score_Intermediate': 1,  # TabPFN score for intermediate classification
        'tabpfn_score_White': -1,        # TabPFN score for white classification (negative)
        #'UMAP2_markers': -0.1,          # UMAP Component 2 based on markers (commented out)
    }

# ==============================================================================
# MAIN EXECUTION PIPELINE
# ==============================================================================

# --- Main Execution ---
def main(protein_data_path = './output/imputed_matrix_2025-11-10.csv', 
         precomputed_lasso_coefs_path='./data/python_input/coefs_from_lasso.pkl', 
         brown_markers_path='./data/python_input/marker.txt', 
         sample_labels_path='./output/meta_data2025-11-10.csv', 
         recompute_lasso_markers=False, 
         output_dir='./data/python_output', compute_shaps=False):
    
    # ==========================================================================
    # SECTION 1: DATA LOADING AND INITIALIZATION
    # ==========================================================================
    
    # --- 0. Load Data ---
    # Load all required datasets including protein expression data, marker genes,
    # sample labels, batch information, and metadata
    protein_df, brown_markers, white_markers, sample_labels, batch_labels, lda_features, markers, tissue, diet= load_example_data(protein_data_path = protein_data_path, precomputed_lasso_coefs_path=precomputed_lasso_coefs_path, 
                      brown_markers_path=brown_markers_path, sample_labels_path=sample_labels_path)

    #print(protein_df.index)
    #print(sample_labels.index)
    # ==========================================================================
    # SECTION 2: FEATURE SELECTION WITH LASSO REGRESSION
    # ==========================================================================
    
    # Alternative to lasso_genes: concatenate brown_markers and white markers
    lasso_genes = brown_markers + white_markers

    # --- Load or Calculate LASSO-Selected Features ---
    # LASSO regression helps identify the most important features for classification
    # by applying L1 regularization to select features with non-zero coefficients
    if recompute_lasso_markers:
        try: 
            with open(f'{output_dir}/lasso.pkl', 'rb') as f:
                lasso_genes = pickle.load(f)
        except FileNotFoundError:
            print("LASSO genes not found, calculating...")
            lasso_genes = lasso_selection(protein_df.join(tissue, how='left'))
            # Save the lasso genes for future use to avoid recalculation
            with open(f'{output_dir}/lasso.pkl', 'wb') as f:
                pickle.dump(lasso_genes, f)


    # --- Generate Base Visualization Plots ---
    # Create comprehensive visualization plots with different feature sets and scaling options
    # These plots provide initial data exploration using PCA, t-SNE, and UMAP
    get_base_plots(protein_df, sample_labels, batch_labels, tissue, diet, 
                   markers=markers, lasso=lasso_genes, all=protein_df.columns.tolist(), scaling=True, output_dir=output_dir)
    get_base_plots(protein_df, sample_labels, batch_labels, tissue, diet, 
                   markers=markers, lasso=lasso_genes, all=protein_df.columns.tolist(), scaling=False, output_dir=output_dir)

    # ==========================================================================
    # SECTION 3: DATA PREPROCESSING AND SETUP
    # ==========================================================================
    
    # --- Filter Available Markers ---
    # Check which marker genes are actually present in the protein expression data
    # This ensures we only work with markers that have corresponding data
    brown_markers_present_in_data = [m for m in brown_markers if m in protein_df.columns]
    white_markers_present_in_data = [m for m in white_markers if m in protein_df.columns]
 
    # --- Initialize Scoring DataFrame ---
    # Create a comprehensive DataFrame to store all scoring results for each sample
    # This will be populated with scores from different machine learning approaches
    all_sample_scores = pd.DataFrame(index=protein_df.index)
    all_sample_scores = all_sample_scores.join(sample_labels, how='left')  # Add sample classification labels
    all_sample_scores = all_sample_scores.join(batch_labels, how='left')   # Add batch information
    # ==========================================================================
    # SECTION 4: PRINCIPAL COMPONENT ANALYSIS (PCA)
    # ==========================================================================
    
    # --- Calculate PCA-based Browning Scores ---
    # PCA provides unsupervised dimensionality reduction to identify major patterns
    # in the protein expression data that may correlate with browning/whitening
    pca_scores, pca_feature_loadings = calculate_pca_browning_score(protein_df, n_components=2)
    all_sample_scores = all_sample_scores.join(pca_scores)
    
    # --- Extract Top Features from PCA Components ---
    # Identify the proteins with highest absolute loadings for each principal component
    # These represent the most influential features for each PC
    pc1_features = pca_feature_loadings['PC1'].abs().sort_values(ascending=False).head(25).index.tolist()
    pc2_features = pca_feature_loadings['PC2'].abs().sort_values(ascending=False).head(25).index.tolist()
    lda_features = pc1_features + pc2_features  # Combined feature set for LDA

    # ==========================================================================
    # SECTION 5: LINEAR DISCRIMINANT ANALYSIS (LDA) MODELING
    # ==========================================================================

    # --- LDA with PCA-selected Features ---
    # Use features identified from PCA components for supervised classification
    # LDA finds linear combinations that best separate brown vs white samples
    print("--- Calculating LDA-based Score (PCA features) ---")
    lda_cv_scores_df, lda_pca_model = calculate_ml_browning_score(
            protein_df.copy()[lda_features], sample_labels, batch_labels,
            positive_class='Brown', negative_class='White',
            clf=LinearDiscriminantAnalysis(), name='lda' 
        )
    all_sample_scores = all_sample_scores.join(lda_cv_scores_df, rsuffix='_pca')

    # --- LDA with LASSO-selected Features ---
    # Use features identified from LASSO regression for supervised classification
    # This approach leverages sparse feature selection for improved interpretability
    print("--- Calculating LDA-based Score (LASSO features) ---")
    lda_cv_scores_df, lda_lasso_model = calculate_ml_browning_score(
            protein_df.copy()[lasso_genes], sample_labels, batch_labels,
            positive_class='Brown', negative_class='White',
            clf=LinearDiscriminantAnalysis(), name='lda' 
        )
    all_sample_scores = all_sample_scores.join(lda_cv_scores_df, rsuffix='_lasso')

    # ==========================================================================
    # SECTION 6: RANDOM FOREST CLASSIFICATION
    # ==========================================================================

    # --- Random Forest with Cross-Validation ---
    # Random Forest provides ensemble-based classification with feature importance
    # Uses all available features and handles non-linear relationships
    print("--- Calculating Random Forest-based Score ---")
    rf_cv_scores_df, rf_model = calculate_ml_browning_score(
            protein_df.copy(), sample_labels, batch_labels,
            positive_class='Brown', negative_class='White',
            clf=RandomForestClassifier(), name='rf'
        )
    all_sample_scores = all_sample_scores.join(rf_cv_scores_df)

    # ==========================================================================
    # SECTION 7: TABPFN CLASSIFICATION
    # ==========================================================================

    # --- TabPFN with Cross-Validation ---
    # TabPFN (Prior-Fitted Networks) is a state-of-the-art method for tabular data
    # that leverages pre-trained neural networks for few-shot learning
    import torch
    print("--- Calculating TabPFN-based Score ---")
    tabpfn_cv_scores_df, tabpfn = calculate_ml_browning_score(
            protein_df.copy(), sample_labels, batch_labels,
            positive_class='Brown', negative_class='White',
            clf=TabPFNClassifier(ignore_pretraining_limits=True, device="cuda" if torch.cuda.is_available() else "cpu"), name='tabpfn'
        )
    all_sample_scores = all_sample_scores.join(tabpfn_cv_scores_df)

    # ==========================================================================
    # SECTION 8: LOGISTIC REGRESSION CLASSIFICATION
    # ==========================================================================

    # --- Logistic Regression with LASSO-selected Features ---
    # Logistic regression provides interpretable linear classification
    # Using LASSO-selected features for sparse and interpretable model
    print("--- Calculating Logistic Regression-based Score (LASSO features) ---")
    lr_cv_scores_df, lr_lasso_model = calculate_ml_browning_score(
            protein_df.copy()[lasso_genes], sample_labels, batch_labels,
            positive_class='Brown', negative_class='White',
            clf=LogisticRegression(max_iter=1000), name='lr'
        )
    all_sample_scores = all_sample_scores.join(lr_cv_scores_df, rsuffix='_lasso')

    # --- Logistic Regression with Marker Genes ---
    # Use biologically-relevant marker genes for classification
    # This approach leverages domain knowledge about browning/whitening markers
    print("--- Calculating Logistic Regression-based Score (Marker genes) ---")
    lr_cv_scores_df, lr_marker_model = calculate_ml_browning_score(
            protein_df.copy()[markers], sample_labels, batch_labels,
            positive_class='Brown', negative_class='White',
            clf=LogisticRegression(max_iter=1000), name='lr'
        )
    all_sample_scores = all_sample_scores.join(lr_cv_scores_df, rsuffix='_markers')

    # ==========================================================================
    # SECTION 9: COMPOSITE SCORE CALCULATION AND FINAL ANALYSIS
    # ==========================================================================

    # --- Calculate Composite Browning Score ---
    # Combine scores from multiple machine learning approaches using weighted average
    # The composite score provides a robust assessment by leveraging multiple models
    composite_score = calculate_composite_score(all_sample_scores.drop(columns=['Status', 'Batch']), weights=None) # if weighting use e.g. example weights
    # --- Save Comprehensive Results ---
    # Store all calculated scores in a CSV file for further analysis and reporting
    all_sample_scores['Composite_Score'] = composite_score
    all_sample_scores.to_csv(f'{output_dir}/all_sample_scores.csv')

    # Display sample results for verification
    print("Sample of calculated scores:")
    print(all_sample_scores.head())

    # --- Correlation Analysis and Visualization ---
    # Analyze relationships between different scoring methods
    # Focus on intermediate samples to understand score consistency
    to_drop = [c for c in all_sample_scores.columns if 'White' in c]  # Remove White-specific scores
    correlation_matrix = all_sample_scores.drop(columns=['Status', 'Batch']+to_drop)[all_sample_scores.Status.eq('Intermediate')].corr()
    
    # Create correlation heatmap for visual analysis
    sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap='coolwarm')
    plt.title("Correlation Matrix of Browning Scores")
    plt.savefig(f"{output_dir}/correlation_matrix.png", dpi=300)
    plt.close()

    # ==========================================================================
    # SECTION 10: MODEL EVALUATION AND FEATURE IMPORTANCE ANALYSIS
    # ==========================================================================

    # --- Feature Importance Analysis for All Models ---
    # Extract and visualize the most important features from each trained model
    # This helps understand which proteins are most predictive of browning
    all_models = [lda_pca_model, lda_lasso_model, rf_model, tabpfn, lr_lasso_model, lr_marker_model]
    model_names = ['LDA PCA', 'LDA LASSO', 'RF', 'TabPFN', 'LR LASSO', 'LR Markers']
    
    
    for model, name in zip(all_models, model_names):
        classifier = model.classifier_
        
        # Extract feature importance based on model type
        if hasattr(classifier, 'feature_importances_'):
            feature_importances = classifier.feature_importances_  # For tree-based models
        elif hasattr(classifier, 'coef_'):
            feature_importances = np.abs(classifier.coef_[0])      # For linear models
        else:
            print(f"Model {name} does not have feature importances or coefficients.")
            continue
        
        # Create feature importance visualization
        feature_names = model.selected_features_
        importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': feature_importances})
        importance_df = importance_df.sort_values(by='Importance', ascending=False).head(25)
        
        plt.figure(figsize=(10, 6))
        sns.barplot(x='Importance', y='Feature', data=importance_df)
        plt.title(f'Top 25 Feature Importances for {name}')
        plt.xlabel('Importance')
        plt.ylabel('Feature')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/feature_importances_{name}.png', dpi=300)
        plt.close()

    # ==========================================================================
    # SECTION 11: SHAP ANALYSIS FOR MODEL INTERPRETABILITY
    # ==========================================================================

    # --- SHAP Analysis for LASSO-selected Features ---
    # SHAP (SHapley Additive exPlanations) provides feature-level explanations
    # for individual predictions, helping understand model decision-making
    print("--- Calculating SHAP values for LASSO-selected features ---")
    if compute_shaps:
        tabpfn_cv_scores_df, tabpfn = calculate_ml_browning_score(
                protein_df.copy()[lasso_genes], sample_labels, batch_labels,
                positive_class='Brown', negative_class='White',
                clf=TabPFNClassifier(ignore_pretraining_limits=True, device='cuda'), name='tabpfn'
            )
        shaps, fit = get_shaps(tabpfn, protein_df.copy()[lasso_genes], sample_labels, file_name="shap_values_lasso", output_dir=output_dir)
        
        # --- SHAP Analysis for Marker Genes ---
        # Analyze SHAP values specifically for biologically-relevant marker genes
        print("--- Calculating SHAP values for marker genes ---")
        tabpfn_cv_scores_df, tabpfn = calculate_ml_browning_score(
                protein_df.copy()[markers], sample_labels, batch_labels,
                positive_class='Brown', negative_class='White',
                clf=TabPFNClassifier(ignore_pretraining_limits=True, device='cuda'), name='tabpfn'
            )
        shaps, fit = get_shaps(tabpfn, protein_df.copy()[markers], sample_labels, file_name="shap_values_markers", output_dir=output_dir)

    # ==========================================================================
    # SECTION 12: FINAL VISUALIZATIONS AND SUMMARY OUTPUTS
    # ==========================================================================

    # --- Generate Score Distribution Plots ---
    # Create comprehensive visualizations showing the distribution of all calculated scores
    # These plots help visualize how well different models separate browning classes
    print("--- Generating score distribution plots ---")
    plot_score_distributions_fix(all_sample_scores.drop(columns=['Status', 'Batch']), all_sample_scores.Status, output_dir=output_dir)

    print("=== AI-BATS Pipeline Completed Successfully ===")
    print(f"Results saved to: script_outputs/")
    print(f"- Comprehensive scores: all_sample_scores.csv")
    print(f"- Visualizations: Various PNG files")
    print(f"- SHAP analysis: SHAP value arrays and plots")
    print(f"- Feature importance: Individual model importance plots")
    print(f"- Correlation analysis: correlation_matrix.png")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="AI-BATS Browning Pipeline")
    parser.add_argument("--protein_data_path", type=str, default='./output/imputed_matrix_2025-11-10.csv')
    parser.add_argument("--precomputed_lasso_coefs_path", type=str, default='./data/python_input/coefs_from_lasso.pkl')
    parser.add_argument("--brown_markers_path", type=str, default='./data/python_input/marker.txt')
    parser.add_argument("--sample_labels_path", type=str, default='./output/meta_data2025-11-10.csv')
    parser.add_argument("--recompute_lasso_markers", type=bool, default=True)
    parser.add_argument("--compute_shaps", type=bool, default=False)
    args = parser.parse_args()
    print(args)
    import time 
    start = time.time()
    main(protein_data_path=args.protein_data_path, precomputed_lasso_coefs_path=args.precomputed_lasso_coefs_path, 
         brown_markers_path=args.brown_markers_path, sample_labels_path=args.sample_labels_path, recompute_lasso_markers=args.recompute_lasso_markers,
         compute_shaps=args.compute_shaps)
    end = time.time()
    print(f"Total time taken: {end - start} seconds")