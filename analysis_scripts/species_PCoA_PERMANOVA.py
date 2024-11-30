import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway
import skbio.stats.ordination as ordination
from skbio.stats.distance import permanova
from skbio import DistanceMatrix
import itertools

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

print(f"Number of species considered: {len(filtered_species_names)}")
print(f"Number of samples: {len(filtered_nests)}")
print("\nNest distribution:")
print(pd.Series(filtered_nests).value_counts())

# Convert to binary presence/absence
abundance_matrix = filtered_species_df.values.astype(float)
abundance_matrix = np.where(abundance_matrix != 0, 1, 0)

# Calculate Bray-Curtis dissimilarity matrix
def bray_curtis_distance(x, y):
    shared = np.sum(np.logical_and(x > 0, y > 0))
    total = np.sum(np.logical_or(x > 0, y > 0))
    if total == 0:
        return 0
    return 1 - (2 * shared) / total

# Calculate distance matrix
dist_matrix = squareform(pdist(abundance_matrix.T, metric=bray_curtis_distance))

# Create sample IDs for the distance matrix
sample_ids = [f"Sample_{i}" for i in range(len(filtered_nests))]

# Create DistanceMatrix object
distance_matrix = DistanceMatrix(dist_matrix, sample_ids)

# Perform PCoA
pcoa = ordination.pcoa(distance_matrix)
pcoa_results = pcoa.samples

# Create a DataFrame with PCoA results and nest information
pcoa_df = pd.DataFrame(data=pcoa_results.values,
                      columns=[f'PC{i+1}' for i in range(pcoa_results.shape[1])])
pcoa_df['Nest'] = filtered_nests

# Calculate variance explained by each axis
variance_explained = pcoa.proportion_explained

# Create PCoA plot
plt.figure(figsize=(10, 8))
colors = sns.color_palette("husl", n_colors=len(np.unique(filtered_nests)))
unique_nests = np.unique(filtered_nests)

# Add sample points
for nest, color in zip(unique_nests, colors):
    mask = pcoa_df['Nest'] == nest
    plt.scatter(pcoa_df.loc[mask, 'PC1'], 
               pcoa_df.loc[mask, 'PC2'],
               c=[color], 
               label=nest,
               alpha=0.7)
    
    # Add centroid as black X
    centroid = pcoa_df.loc[mask, ['PC1', 'PC2']].mean()
    plt.scatter(centroid['PC1'], centroid['PC2'], 
               c='black', marker='x', s=100, 
               linewidth=2)

plt.xlabel(f'PCoA 1 ({variance_explained[0]:.1%} explained)')
plt.ylabel(f'PCoA 2 ({variance_explained[1]:.1%} explained)')
plt.title('PCoA of Prey Species Composition in Guma\n(Bray-Curtis dissimilarity)')
plt.legend(title='Nest', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()

# Save the PCoA plot
save_dir = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics'
plt.savefig(f'{save_dir}/pcoa_guma.png', dpi=300, bbox_inches='tight')
plt.show()

# Perform simple PERMANOVA using custom function
def custom_permanova(distance_matrix, groups, n_perm=999):
    """
    Simple PERMANOVA implementation
    """
    # Convert to numpy array if not already
    distance_matrix = np.array(distance_matrix)
    
    # Calculate actual F-statistic
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)
    
    # Total sum of squares
    total_ss = np.sum(distance_matrix ** 2) / len(groups)
    
    # Within-group sum of squares
    within_ss = 0
    for group in unique_groups:
        mask = groups == group
        if np.sum(mask) > 1:  # Only if more than one sample in group
            group_distances = distance_matrix[mask][:, mask]
            within_ss += np.sum(group_distances ** 2) / np.sum(mask)
    
    # Between-group sum of squares
    between_ss = total_ss - within_ss
    
    # Degrees of freedom
    df_between = n_groups - 1
    df_within = len(groups) - n_groups
    
    # F-statistic
    f_stat = (between_ss / df_between) / (within_ss / df_within)
    
    # Permutation test
    perm_f_stats = np.zeros(n_perm)
    for i in range(n_perm):
        perm_groups = np.random.permutation(groups)
        
        # Calculate within-group sum of squares for permuted groups
        perm_within_ss = 0
        for group in unique_groups:
            mask = perm_groups == group
            if np.sum(mask) > 1:
                group_distances = distance_matrix[mask][:, mask]
                perm_within_ss += np.sum(group_distances ** 2) / np.sum(mask)
        
        # Calculate permuted F-statistic
        perm_between_ss = total_ss - perm_within_ss
        perm_f_stats[i] = (perm_between_ss / df_between) / (perm_within_ss / df_within)
    
    # Calculate p-value
    p_value = np.mean(perm_f_stats >= f_stat)
    
    return {'F_statistic': f_stat, 'p_value': p_value, 'permutations': n_perm}

# Perform PERMANOVA
permanova_results = custom_permanova(dist_matrix, filtered_nests, n_perm=999)

print("\nPERMANOVA Results:")
print(f"F-statistic: {permanova_results['F_statistic']:.4f}")
print(f"p-value: {permanova_results['p_value']:.4f}")
print(f"Number of permutations: {permanova_results['permutations']}")

# Pairwise PERMANOVA
pairwise_results = []

for nest1, nest2 in itertools.combinations(unique_nests, 2):
    # Get indices for the two nests
    mask1 = filtered_nests == nest1
    mask2 = filtered_nests == nest2
    combined_mask = mask1 | mask2
    
    # Extract relevant distances and grouping
    sub_dist = dist_matrix[combined_mask][:, combined_mask]
    sub_groups = filtered_nests[combined_mask]
    
    # Perform PERMANOVA
    result = custom_permanova(sub_dist, sub_groups, n_perm=999)
    
    pairwise_results.append({
        'Nest1': nest1,
        'Nest2': nest2,
        'F-statistic': result['F_statistic'],
        'p-value': result['p_value']
    })

# Create DataFrame of pairwise results
pairwise_df = pd.DataFrame(pairwise_results)

print("\nPairwise PERMANOVA Results:")
print(pairwise_df.to_string(index=False))