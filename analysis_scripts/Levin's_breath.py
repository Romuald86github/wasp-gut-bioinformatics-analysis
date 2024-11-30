import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Load the presence-absence data
data_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/presence_absence_IK_data.xlsx'
data_df = pd.read_excel(data_path)

# Extract metadata (Nest and Location) and species data
metadata_df = data_df.iloc[0:2, 1:]  # First two rows contain Nest and Location
species_df = data_df.iloc[2:, 1:]    # Species data starts from row 3
species_names = data_df.iloc[2:, 0]   # Species names in first column

# Convert metadata to arrays for filtering
nests = metadata_df.iloc[0].values
locations = metadata_df.iloc[1].values

# Create mask for Guma location and excluding nests S and T
location_mask = (locations == 'Guma')
nest_mask = ~np.isin(nests, ['S', 'T'])
combined_mask = location_mask & nest_mask

# Apply the mask to metadata and species data
filtered_nests = nests[combined_mask]
filtered_species_df = species_df.iloc[:, combined_mask]

# Function to check if a name is a species (has genus and species name)
def is_species_name(name):
    parts = name.split()
    if len(parts) == 2:
        return parts[0][0].isupper() and parts[1].islower()
    return False

# Filter to keep only true species names
species_mask = species_names.apply(is_species_name)
filtered_species_df = filtered_species_df[species_mask]
filtered_species_names = species_names[species_mask]

# Convert to binary presence/absence
abundance_matrix = filtered_species_df.values.astype(float)
abundance_matrix = np.where(abundance_matrix != 0, 1, 0)

# Calculate Levins' niche breadth for each nest
def levins_breadth(sample):
    n_used = np.sum(sample)  # number of prey species used
    if n_used == 0:
        return 0
    # For presence/absence data, all used resources have equal proportions
    p = 1.0 / n_used
    B = 1 / (n_used * p * p)  # Levins' B
    # Standardized Levins' B (BA)
    n_available = len(sample)  # total number of prey species available
    BA = (B - 1) / (n_available - 1) if n_available > 1 else 0
    return BA

# Calculate Levins' breadth for each sample
levins_indices = [levins_breadth(abundance_matrix[:, i]) for i in range(abundance_matrix.shape[1])]

# Create DataFrame with nest information and Levins' indices
niche_df = pd.DataFrame({
    'Nest': filtered_nests,
    'Levins_B': levins_indices
})

# Calculate summary statistics
summary_stats = pd.DataFrame({
    'Statistic': ['Mean', 'Median', 'Std Dev', 'Min', 'Max'],
    'Value': [
        np.mean(levins_indices),
        np.median(levins_indices),
        np.std(levins_indices),
        np.min(levins_indices),
        np.max(levins_indices)
    ]
})

print("Summary Statistics of Standardized Levins' Niche Breadth:")
print(summary_stats)
print("\nSample size:", len(levins_indices))

# Create histogram of Levins' breadth
plt.figure(figsize=(10, 6))
sns.histplot(data=niche_df, x='Levins_B', bins=20)
plt.axvline(np.mean(levins_indices), color='red', linestyle='--', label='Mean')
plt.axvline(np.median(levins_indices), color='green', linestyle='--', label='Median')
plt.title("Distribution of Levins' Niche Breadth in Guma\n(excluding nests S and T)")
plt.xlabel("Standardized Levins' Niche Breadth (BA)")
plt.ylabel('Count')
plt.legend()
plt.tight_layout()

# Save the histogram
save_dir = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics'
plt.savefig(f'{save_dir}/levins_breadth_distribution_guma.png', dpi=300, bbox_inches='tight')
plt.show()

