import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

# Calculate Shannon diversity for each sample
def shannon_diversity(sample):
    # Count number of present species
    n_present = np.sum(sample)
    if n_present == 0:
        return 0
    
    # For presence/absence data, all present species have equal relative abundance
    p = 1.0 / n_present
    
    # Calculate Shannon index
    # H = -sum(p_i * ln(p_i)) where all p_i are equal
    return -n_present * (p * np.log(p))

# Calculate Shannon diversity for each sample
shannon_indices = [shannon_diversity(abundance_matrix[:, i]) for i in range(abundance_matrix.shape[1])]

# Create a dataframe with nest information and Shannon indices
diversity_df = pd.DataFrame({
    'Nest': filtered_nests,
    'Shannon_Index': shannon_indices
})

# Calculate summary statistics
summary_stats = pd.DataFrame({
    'Statistic': ['Mean', 'Median', 'Std Dev', 'Min', 'Max'],
    'Value': [
        np.mean(shannon_indices),
        np.median(shannon_indices),
        np.std(shannon_indices),
        np.min(shannon_indices),
        np.max(shannon_indices)
    ]
})

print("Summary Statistics of Shannon Diversity:")
print(summary_stats)
print("\nSample size:", len(shannon_indices))

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), height_ratios=[2, 1])

# Histogram
sns.histplot(data=diversity_df, x='Shannon_Index', bins=20, ax=ax1)
ax1.axvline(np.mean(shannon_indices), color='red', linestyle='--', label='Mean')
ax1.axvline(np.median(shannon_indices), color='green', linestyle='--', label='Median')
ax1.set_title('Distribution of Shannon Diversity Indices in Guma\n(excluding nests S and T)')
ax1.set_xlabel('Shannon Diversity Index')
ax1.set_ylabel('Count')
ax1.legend()

# Boxplot
sns.boxplot(data=diversity_df, x='Shannon_Index', ax=ax2)
ax2.set_title('Boxplot of Shannon Diversity Indices')
ax2.set_xlabel('Shannon Diversity Index')

# Add individual points over the boxplot
sns.stripplot(data=diversity_df, x='Shannon_Index', color='red', size=4, alpha=0.3, ax=ax2)

# Adjust layout
plt.tight_layout()

# Save the plot
save_dir = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics'
plt.savefig(f'{save_dir}/shannon_diversity_distribution_guma.png', dpi=300, bbox_inches='tight')
plt.show()

# Create boxplot by nest
plt.figure(figsize=(12, 6))
sns.boxplot(data=diversity_df, x='Nest', y='Shannon_Index')
sns.stripplot(data=diversity_df, x='Nest', y='Shannon_Index', color='red', size=4, alpha=0.3)
plt.title('Shannon Diversity Indices by Nest in Guma\n(excluding nests S and T)')
plt.xlabel('Nest')
plt.ylabel('Shannon Diversity Index')
plt.xticks(rotation=45)
plt.tight_layout()

# Save the nest-wise plot
plt.savefig(f'{save_dir}/shannon_diversity_by_nest_guma.png', dpi=300, bbox_inches='tight')
plt.show()