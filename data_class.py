import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from typing import Tuple, List

class DataLoader:
    """
    A class for loading and processing protein data with configurable file paths.
    """
    
    def __init__(self, config_path: str = "config.yaml"):
        """
        Initialize the DataLoader with configuration from YAML file.
        
        Args:
            config_path: Path to the YAML configuration file
        """
        self.config_path = Path(config_path)
        self.config = self._load_config()
        
        # Store loaded data
        self.protein_df = None
        self.brown_markers_list = None
        self.white_markers_list = None
        self.sample_labels_series = None
        self.batch_labels_series = None
        self.lda_features = None
        self.brown_markers_from_file = None
        self.train_idx = []
        self.test_idx = []
    
    def _load_config(self) -> dict:
        """Load configuration from YAML file."""
        try:
            with open(self.config_path, 'r') as file:
                return yaml.safe_load(file)
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML configuration: {e}")
    
    def load_data(self) -> Tuple[pd.DataFrame, List, List, pd.Series, pd.Series, List, List, List, List]:
        """
        Load and process all data according to configuration.
        
        Returns:
            Tuple containing: protein_df, brown_markers_list, white_markers_list, 
            sample_labels_series, batch_labels_series, lda_features, 
            brown_markers_from_file, train_idx, test_idx
        """
        # Load protein data
        self.protein_df = pd.read_csv(
            self.config['data_paths']['protein_data'], 
            index_col=0, 
            header=0
        ).T
        
        # Load coefficients
        coefs = pd.read_pickle(self.config['data_paths']['coefficients']).T
        
        # Load brown markers from file
        self.brown_markers_from_file = pd.read_csv(
            self.config['data_paths']['markers'], 
            sep="\t", 
            index_col=0, 
            header=0
        )["GeneName"].to_list()
        self.brown_markers_from_file = [
            i for i in self.brown_markers_from_file 
            if i in self.protein_df.columns
        ]
        
        # Load metadata
        metadata = pd.read_csv(
            self.config['data_paths']['metadata'], 
            index_col=0, 
            header=0
        )
        
        # Process batch information
        batchd = {s: s.split('-')[0] for s in metadata.file_name}
        batch = [batchd[s] for s in self.protein_df.index]
        
        # Extract markers based on coefficient threshold
        threshold = self.config['processing']['coefficient_threshold']
        wc = [i for j, i in zip(coefs.eWAT, coefs.index) if abs(j) > threshold]
        bc = [i for j, i in zip(coefs.BAT, coefs.index) if abs(j) > threshold]
        
        self.brown_markers_list = self.brown_markers_from_file + bc
        self.white_markers_list = wc
        
        # Process sample labels
        sample_labels = pd.read_csv(
            self.config['data_paths']['metadata'], 
            index_col=0, 
            header=0
        )
        
        tissue_mapping = self.config['tissue_mapping']
        sl = [tissue_mapping[i] for i in sample_labels.tissue.to_list()]
        
        self.sample_labels_series = pd.Series(
            sl, 
            index=sample_labels.index, 
            name="Status"
        )
        self.batch_labels_series = pd.Series(
            batch, 
            index=sample_labels.index, 
            name="Batch"
        )
        
        # Combine features for LDA
        self.lda_features = self.brown_markers_list + self.white_markers_list
        
        return (
            self.protein_df, 
            self.brown_markers_list, 
            self.white_markers_list, 
            self.sample_labels_series, 
            self.batch_labels_series, 
            self.lda_features, 
            self.brown_markers_from_file, 
            self.train_idx, 
            self.test_idx
        )
    
    def get_protein_data(self) -> pd.DataFrame:
        """Get the loaded protein dataframe."""
        if self.protein_df is None:
            self.load_data()
        return self.protein_df
    
    def get_markers(self) -> Tuple[List, List]:
        """Get brown and white marker lists."""
        if self.brown_markers_list is None or self.white_markers_list is None:
            self.load_data()
        return self.brown_markers_list, self.white_markers_list
    
    def get_sample_labels(self) -> pd.Series:
        """Get sample labels series."""
        if self.sample_labels_series is None:
            self.load_data()
        return self.sample_labels_series
    
    def get_batch_labels(self) -> pd.Series:
        """Get batch labels series."""
        if self.batch_labels_series is None:
            self.load_data()
        return self.batch_labels_series


def load_example_data():
    """
    Legacy function for backward compatibility.
    Creates a DataLoader instance and returns the same data structure.
    """
    loader = DataLoader()
    return loader.load_data()