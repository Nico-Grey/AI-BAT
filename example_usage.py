# Example usage of the DataLoader class
from data_functions import DataLoader

# Create a DataLoader instance with default config
loader = DataLoader()

# Load all data
data = loader.load_data()
protein_df, brown_markers, white_markers, sample_labels, batch_labels, lda_features, brown_markers_from_file, train_idx, test_idx = data

# Or use individual methods
protein_data = loader.get_protein_data()
brown_markers, white_markers = loader.get_markers()
sample_labels = loader.get_sample_labels()
batch_labels = loader.get_batch_labels()

# For backward compatibility, the original function still works
from data_functions import load_example_data
legacy_data = load_example_data()

print("DataLoader class created successfully!")
print(f"Protein data shape: {protein_data.shape}")
print(f"Number of brown markers: {len(brown_markers)}")
print(f"Number of white markers: {len(white_markers)}")
