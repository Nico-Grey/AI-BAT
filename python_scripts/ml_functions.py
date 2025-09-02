import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from tabpfn import TabPFNClassifier
# from tabpfn import TabPFNClassifier
# from analysis_functions import get_tabpfn_score

def lasso_selection(protein_df):
    """Performs LASSO regression to select features based on their coefficients."""

    X = protein_df.drop(columns=['Tissue'])
    y = LabelEncoder().fit_transform(protein_df['Tissue'])  


    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    lasso = LogisticRegression(penalty='l1', solver='saga', max_iter=25000)
    lasso.fit(X_scaled, y)

    selected_features = pd.DataFrame(lasso.coef_.T, index=X.columns)
    selected_features = selected_features.loc[(selected_features != 0).any(axis=1)]

    print(f"Selected {selected_features.shape} features using LASSO.")
    
    return selected_features.index.tolist()

"""def calculate_ml_browning_score(protein_df, sample_labels_full, batch_labels_full,
                                 positive_class='Brown', negative_class='White',feature_selection_pca_loadings=None,
                                 clf=LinearDiscriminantAnalysis(), name = 'lda'):
    protein_expr = protein_df.join(sample_labels_full, how='inner')
    train = protein_expr[protein_expr['Status'].isin([positive_class, negative_class])]
    #test = protein_expr[~protein_expr['Status'].isin([positive_class, negative_class])]
    if feature_selection_pca_loadings is not None and name == 'lda':
        # Filter features based on PCA loadings
        selected_features = [f for f in feature_selection_pca_loadings if f in train.columns]
    else:
        selected_features = train.columns[:-1]  # Exclude the 'Status' column
    scaler = StandardScaler()
    clf.fit(scaler.fit_transform(train[selected_features]), train['Status'])
    lda_scores = pd.DataFrame(clf.predict_proba(scaler.transform(protein_df[selected_features])),
                              index=protein_df.index,
                              columns=[f'{name}_browning_score_{cls}' for cls in clf.classes_])
    return lda_scores, clf"""

from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, TransformerMixin


class MLBrowningPipeline:
    """Pipeline wrapper that tracks input features and provides browning scores."""
    
    def __init__(self, positive_class='Brown', negative_class='White', 
                 feature_selection_pca_loadings=None, clf=None, name='lda'):
        self.positive_class = positive_class
        self.negative_class = negative_class
        self.name = name
        
        # Create pipeline
        self.pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', clf if clf is not None else LinearDiscriminantAnalysis())
        ])
        
        self.is_fitted_ = False
        self.selected_features = None
    
    def fit(self, protein_df, sample_labels_full):
        """Fit the pipeline on training data."""

        protein_expr = protein_df.join(sample_labels_full, how='inner')
        train = protein_expr[protein_expr['Status'].isin([self.positive_class, self.negative_class])]
        
        X_train = train.drop(columns=['Status'])
        y_train = train['Status']
        self.selected_features = X_train.columns.tolist()

        self.pipeline.fit(X_train, y_train)
        self.is_fitted_ = True

        return self
    
    def predict_proba(self, protein_df):
        """Get probability predictions for all samples."""
        if not self.is_fitted_:
            raise ValueError("Pipeline must be fitted before prediction")
        
        # Only use features that were selected during training
        X = protein_df#[self.selected_features_]
        probas = self.pipeline.predict_proba(X)
        
        scores = pd.DataFrame(
            probas,
            index=protein_df.index,
            columns=[f'{self.name}_browning_score_{cls}' for cls in self.pipeline.named_steps['classifier'].classes_]
        )
        return scores
    
    @property
    def selected_features_(self):
        """Get the list of selected features."""
        if not self.is_fitted_:
            return None
        return self.selected_features
    
    @property
    def classifier_(self):
        """Get the fitted classifier."""
        if not self.is_fitted_:
            return None
        return self.pipeline.named_steps['classifier']

def calculate_ml_browning_score(protein_df, sample_labels_full, batch_labels_full,
                                positive_class='Brown', negative_class='White',
                                feature_selection_pca_loadings=None,
                                clf=LinearDiscriminantAnalysis(), name='lda'):
    """
    Calculate ML browning scores using a pipeline approach.
    Returns scores and a pipeline object that tracks input features.
    """
    # Create and fit the pipeline
    ml_pipeline = MLBrowningPipeline(
        positive_class=positive_class,
        negative_class=negative_class,
        feature_selection_pca_loadings=feature_selection_pca_loadings,
        clf=clf,
        name=name
    )
    
    ml_pipeline.fit(protein_df, sample_labels_full)
    scores = ml_pipeline.predict_proba(protein_df)
    
    return scores, ml_pipeline

def get_tabpfn_score(protein_df, sample_labels):
    cv = LeaveOneOut().split(protein_df)
    scores = pd.DataFrame(index=protein_df.index,columns=[f'tabpfn_score_{class_name}' for class_name in sample_labels.unique()])
    clf = TabPFNClassifier(ignore_pretraining_limits=True, device='cuda')  # Ensure you have the correct device set up
    for i, (train_index, test_index) in enumerate(cv):
        train_data = protein_df.iloc[train_index]
        test_data = protein_df.iloc[test_index]
        t_index = protein_df.index[test_index]
        principal_components = clf.fit(train_data, sample_labels.iloc[train_index])
        test_pc = clf.predict_proba(test_data)
        for i, class_name in enumerate(sample_labels.unique()):
            scores.loc[t_index, f'tabpfn_score_{class_name}' ] = test_pc[:, i]

    return scores, clf


def get_shaps(tabpfn, protein_df, sample_labels, file_name="shap_values"):
    """
    Get SHAP values for the TabPFN model.
    """
    from tabpfn_extensions import interpretability
    feature_names = protein_df.columns.tolist()
    X_test = protein_df.values
    shap_values = interpretability.shap.get_shap_values(
    estimator=tabpfn,
    test_x=X_test,
    attribute_names=feature_names,
    algorithm="permutation",
    max_evals = 1500,)

    np.save(f"script_outputs/{file_name}.npy", shap_values)

# Create visualization
    fig = interpretability.shap.plot_shap(shap_values)
    
    return shap_values, fig
