# browning_pipeline.py - Advanced Machine Learning Pipeline for Adipose Tissue Analysis


NOTE - This is not up to date :::TODO:::

## Overview

browning_pipeline.py is a comprehensive Python-based machine learning pipeline designed for analyzing protein expression data to classify and score brown versus white adipose tissue. The pipeline implements multiple state-of-the-art machine learning algorithms, feature selection methods, and interpretability tools to provide robust browning scores for biological samples.

## Features

### 🤖 Machine Learning Models
- **Linear Discriminant Analysis (LDA)**: Classical linear classification
- **Logistic Regression**: With LASSO regularization
- **Random Forest**: Ensemble method with feature importance
- **TabPFN**: Neural network-based tabular predictor
- **Cross-validation**: Stratified k-fold validation for robust scoring

### 🎯 Feature Selection
- **LASSO Regularization**: Automatic feature selection
- **PCA-based Selection**: Top principal component loadings
- **Marker-based Analysis**: Predefined brown/white markers
- **Feature Importance Ranking**: Model-specific importance scores

### 📊 Interpretability & Visualization
- **SHAP Analysis**: Model-agnostic explanations
- **Feature Importance Plots**: Visual ranking of predictive features
- **Score Distributions**: Brown vs. white vs. intermediate comparisons
- **Correlation Analysis**: Score relationship visualization

### 🔬 Dimensionality Reduction
- **PCA**: Principal Component Analysis
- **t-SNE**: t-distributed Stochastic Neighbor Embedding
- **UMAP**: Uniform Manifold Approximation and Projection

## Dependencies

### Core Libraries
```python
# Machine Learning
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, cross_val_score

# Deep Learning
from tabpfn import TabPFNClassifier

# Data Processing
import numpy as np
import pandas as pd

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

# Dimensionality Reduction
import umap

# Interpretability
import shap

# Statistics
from scipy.stats import mannwhitneyu

# Utilities
from pathlib import Path
import pickle
```

## Installation

### Prerequisites
```bash
# Create conda environment
conda create -n aibats python=3.8
conda activate aibats

# Install core packages
pip install scikit-learn pandas numpy matplotlib seaborn

# Install specialized packages
pip install umap-learn tabpfn shap

# For GPU support (optional)
pip install torch torchvision  # for TabPFN GPU acceleration
```

### Custom Modules
The pipeline requires custom modules located in the same directory:
- `data_functions.py`: Data loading and preprocessing
- `ml_functions.py`: Machine learning utilities
- `scoring_functions.py`: Scoring and evaluation functions
- `visualisation_functions.py`: Plotting and visualization
- `analysis_functions.py`: Statistical analysis tools

## Usage

### Basic Execution
```python
python browning_pipeline.py
```

### Programmatic Usage
```python
from browning_pipeline import *

# Load your data
protein_df, brown_markers, white_markers, sample_labels, batch_labels = load_your_data()

# Run specific analysis
pca_scores = calculate_pca_browning_score(protein_df, n_components=2)
ml_scores = calculate_ml_browning_score(protein_df, sample_labels, batch_labels)
```

## Data Format Requirements

### Input Data Structure

#### Protein Expression DataFrame
```python
# protein_df: pandas DataFrame
# - Rows: Samples
# - Columns: Protein/gene identifiers
# - Values: Expression levels or normalized intensities
```

#### Sample Labels
```python
# sample_labels: pandas Series
# - Index: Sample identifiers (matching protein_df index)
# - Values: 'Brown', 'White', or 'Intermediate'
```

#### Batch Labels
```python
# batch_labels: pandas Series  
# - Index: Sample identifiers
# - Values: Batch/experiment identifiers
```

#### Marker Lists
```python
# brown_markers: list of strings
# - Protein/gene identifiers associated with brown adipose tissue

# white_markers: list of strings  
# - Protein/gene identifiers associated with white adipose tissue
```

### Example Data Structure
```python
import pandas as pd

# Protein expression data
protein_df = pd.DataFrame({
    'UCP1': [8.5, 2.1, 7.8, 1.9],
    'CIDEA': [6.2, 1.5, 5.9, 1.2],
    'COX8B': [4.3, 0.8, 4.1, 0.9],
    # ... more proteins
}, index=['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4'])

# Sample classifications
sample_labels = pd.Series(['Brown', 'White', 'Brown', 'White'], 
                         index=protein_df.index, name='Status')

# Batch information
batch_labels = pd.Series(['Batch_1', 'Batch_1', 'Batch_2', 'Batch_2'],
                        index=protein_df.index, name='Batch')
```

## Analysis Pipeline

### 1. Data Loading and Preprocessing
```python
# Load example data or your custom data
protein_df, brown_markers, white_markers, sample_labels, batch_labels, \
lda_features, markers, tissue, diet = load_example_data()
```

### 2. Feature Selection (LASSO)
```python
# Perform LASSO feature selection
try:
    with open('script_outputs/lasso.pkl', 'rb') as f:
        lasso_genes = pickle.load(f)
except FileNotFoundError:
    lasso_genes = lasso_selection(protein_df.join(tissue, how='left'))
    with open('script_outputs/lasso.pkl', 'wb') as f:
        pickle.dump(lasso_genes, f)
```

### 3. PCA-Based Scoring
```python
# Calculate PCA scores
pca_scores, pca_feature_loadings = calculate_pca_browning_score(
    protein_df, n_components=2
)

# Extract top features from PC loadings
pc1_features = pca_feature_loadings['PC1'].abs().sort_values(ascending=False).head(25).index.tolist()
pc2_features = pca_feature_loadings['PC2'].abs().sort_values(ascending=False).head(25).index.tolist()
```

### 4. Machine Learning Models

#### Linear Discriminant Analysis
```python
# LDA with PCA features
lda_pca_scores, lda_pca_model = calculate_ml_browning_score(
    protein_df[lda_features], sample_labels, batch_labels,
    positive_class='Brown', negative_class='White',
    clf=LinearDiscriminantAnalysis(), name='lda'
)

# LDA with LASSO features
lda_lasso_scores, lda_lasso_model = calculate_ml_browning_score(
    protein_df[lasso_genes], sample_labels, batch_labels,
    positive_class='Brown', negative_class='White',
    clf=LinearDiscriminantAnalysis(), name='lda'
)
```

#### Random Forest
```python
rf_scores, rf_model = calculate_ml_browning_score(
    protein_df, sample_labels, batch_labels,
    positive_class='Brown', negative_class='White',
    clf=RandomForestClassifier(), name='rf'
)
```

#### TabPFN (Neural Network)
```python
tabpfn_scores, tabpfn_model = calculate_ml_browning_score(
    protein_df, sample_labels, batch_labels,
    positive_class='Brown', negative_class='White',
    clf=TabPFNClassifier(ignore_pretraining_limits=True, device='cuda'),
    name='tabpfn'
)
```

#### Logistic Regression
```python
# With LASSO features
lr_lasso_scores, lr_lasso_model = calculate_ml_browning_score(
    protein_df[lasso_genes], sample_labels, batch_labels,
    positive_class='Brown', negative_class='White',
    clf=LogisticRegression(max_iter=1000), name='lr'
)

# With marker features
lr_marker_scores, lr_marker_model = calculate_ml_browning_score(
    protein_df[markers], sample_labels, batch_labels,
    positive_class='Brown', negative_class='White',
    clf=LogisticRegression(max_iter=1000), name='lr'
)
```

### 5. Composite Scoring
```python
# Define weights for different scoring methods
example_weights = {
    'PC2_markers': 0.1,
    'PC2': -0.1,
    'LogisticRegression_probability_score': 0.4,
    'RandomForest_probability_score': 0.5,
    'tabpfn_score_Brown': 3,
    'tabpfn_score_Intermediate': 1,
    'tabpfn_score_White': -1,
}

# Calculate composite score
composite_score = calculate_composite_score(
    all_sample_scores.drop(columns=['Status', 'Batch']), 
    example_weights
)
```

## Output Files

### Scoring Results
- `all_sample_scores.csv`: Complete scoring matrix for all samples
- Individual model scores and probability estimates
- Composite browning scores

### Visualizations
- `correlation_matrix.png`: Score correlation heatmap
- `feature_importances_[MODEL].png`: Feature importance plots for each model
- `[score]_distribution.png`: Score distribution plots by tissue type
- `Composite_Score_distribution.png`: Final composite score distributions

### Model Artifacts
- `lasso.pkl`: LASSO-selected features
- `shap_values_[features].npy`: SHAP explanation values
- Model objects with fitted parameters

### SHAP Analysis
- `Figure_1_shap_lasso.png`: SHAP summary plot for LASSO features
- `Figure_2_shap_lasso.png`: SHAP dependence plots
- `Figure_2_shap_marker.png`: SHAP analysis for marker features

## Configuration

### Model Parameters

#### LASSO Feature Selection
```python
lasso_params = {
    'alpha': 0.01,  # Regularization strength
    'max_iter': 1000,
    'selection': 'cyclic'
}
```

#### Random Forest
```python
rf_params = {
    'n_estimators': 100,
    'max_depth': None,
    'min_samples_split': 2,
    'random_state': 42
}
```

#### TabPFN
```python
tabpfn_params = {
    'ignore_pretraining_limits': True,
    'device': 'cuda',  # or 'cpu'
    'N_ensemble_configurations': 32
}
```

### Cross-Validation Settings
```python
cv_params = {
    'cv': 5,  # 5-fold cross-validation
    'stratify': True,  # Maintain class balance
    'random_state': 42
}
```

## Advanced Features

### SHAP Interpretability
```python
# Generate SHAP explanations
shaps, fit = get_shaps(
    tabpfn_model, 
    protein_df[lasso_genes], 
    sample_labels, 
    file_name="shap_values_lasso"
)
```

### Custom Scoring Weights
```python
# Define custom composite score weights
custom_weights = {
    'model1_score': 0.3,
    'model2_score': 0.4,
    'model3_score': 0.3
}

composite_score = calculate_composite_score(scores_df, custom_weights)
```

### Base Plots Generation
```python
# Generate dimensionality reduction plots
get_base_plots(
    protein_df, sample_labels, batch_labels, tissue, diet,
    markers=markers, lasso=lasso_genes, all=protein_df.columns.tolist(),
    scaling=True
)
```

## Performance Optimization

### GPU Acceleration
```python
# For TabPFN with CUDA support
tabpfn = TabPFNClassifier(
    ignore_pretraining_limits=True,
    device='cuda'  # Requires CUDA-compatible GPU
)
```

### Memory Management
```python
# For large datasets
import gc

# Clear memory after intensive operations
del large_dataframe
gc.collect()

# Use chunking for very large datasets
def process_in_chunks(df, chunk_size=1000):
    for chunk in pd.read_csv('large_file.csv', chunksize=chunk_size):
        yield process_chunk(chunk)
```

### Parallel Processing
```python
from joblib import Parallel, delayed

# Parallel cross-validation
scores = Parallel(n_jobs=-1)(
    delayed(cross_val_score)(model, X, y, cv=cv) 
    for model in models
)
```

## Troubleshooting

### Common Issues

#### CUDA/GPU Problems
```python
# Check CUDA availability
import torch
print(f"CUDA available: {torch.cuda.is_available()}")

# Fallback to CPU
tabpfn = TabPFNClassifier(device='cpu')
```

#### Memory Errors
```python
# Reduce dataset size
sampled_df = protein_df.sample(frac=0.8, random_state=42)

# Use float32 instead of float64
protein_df = protein_df.astype('float32')
```

#### Missing Features
```python
# Handle missing markers
available_markers = [m for m in markers if m in protein_df.columns]
if len(available_markers) < len(markers):
    print(f"Warning: {len(markers) - len(available_markers)} markers not found")
```

### Debugging
```python
# Enable verbose output
import logging
logging.basicConfig(level=logging.DEBUG)

# Check data shapes
print(f"Protein data shape: {protein_df.shape}")
print(f"Sample labels shape: {sample_labels.shape}")
print(f"Batch labels shape: {batch_labels.shape}")
```

## Model Validation

### Cross-Validation Metrics
- **Accuracy**: Overall classification accuracy
- **AUC-ROC**: Area under the ROC curve
- **Precision/Recall**: Class-specific performance
- **F1-Score**: Harmonic mean of precision and recall

### Statistical Testing
```python
# Mann-Whitney U test for score differences
from scipy.stats import mannwhitneyu

brown_scores = scores[sample_labels == 'Brown']
white_scores = scores[sample_labels == 'White']
statistic, p_value = mannwhitneyu(brown_scores, white_scores)
```

## Integration with R Pipeline

The Python pipeline integrates seamlessly with the R-based BATpipeline.R:

```r
# From R, call Python pipeline
system("python script_modularized/browning_pipeline.py")

# Or using reticulate
library(reticulate)
source_python("script_modularized/browning_pipeline.py")
```

## Best Practices

### Data Quality
1. **Check for missing values**: Handle NaN/NA appropriately
2. **Validate sample matching**: Ensure consistent sample IDs across datasets
3. **Batch effect assessment**: Visualize batch effects before modeling

### Model Selection
1. **Feature selection**: Use LASSO or other methods to avoid overfitting
2. **Cross-validation**: Always use proper CV for model evaluation
3. **Ensemble methods**: Combine multiple models for robust predictions

### Interpretability
1. **SHAP analysis**: Always generate explanations for model predictions
2. **Feature importance**: Understand which proteins drive classifications
3. **Biological validation**: Ensure results align with biological knowledge

## Support and Maintenance

### Version Compatibility
- Python: 3.7+
- scikit-learn: 0.24+
- pandas: 1.3+
- NumPy: 1.20+

### Performance Monitoring
```python
import time
import psutil

# Monitor execution time
start_time = time.time()
# ... your code ...
execution_time = time.time() - start_time

# Monitor memory usage
memory_usage = psutil.Process().memory_info().rss / 1024 / 1024  # MB
```
