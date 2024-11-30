import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
import skbio.stats.ordination as ordination
from skbio import DistanceMatrix

# Load and process data
data_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/presence_absence_IK_data.xlsx'
data_df = pd.read_excel(data_path)

metadata_df = data_df.iloc[0:2, 1:]
species_df = data_df.iloc[2:, 1:]
original_names = data_df.iloc[2:, 0]

# Filter for Guma and exclude nests S and T
nests = metadata_df.iloc[0].values
locations = metadata_df.iloc[1].values
location_mask = (locations == 'Guma')
nest_mask = ~np.isin(nests, ['S', 'T'])
combined_mask = location_mask & nest_mask

filtered_nests = nests[combined_mask]
filtered_species_df = species_df.iloc[:, combined_mask]

# Extract genera
def get_genus(name):
    return name.split()[0]

genera_names = original_names.apply(get_genus)

# Create genera presence/absence matrix
genera_to_rows = {}
for idx, genus in enumerate(genera_names):
    if genus not in genera_to_rows:
        genera_to_rows[genus] = []
    genera_to_rows[genus].append(idx)

# Convert to binary presence/absence
abundance_matrix = filtered_species_df.values.astype(float)
abundance_matrix = np.where(abundance_matrix != 0, 1, 0)

# Create genera presence matrix
unique_genera = sorted(genera_to_rows.keys())
genera_matrix = np.zeros((len(unique_genera), abundance_matrix.shape[1]))
for i, genus in enumerate(unique_genera):
    row_indices = genera_to_rows[genus]
    genera_matrix[i, :] = np.any(abundance_matrix[row_indices, :], axis=0)

# 1. Plot genera accumulation curve
def accumulate_genera(genera_matrix, n_permutations=100):
    n_samples = genera_matrix.shape[1]
    accumulation_curves = np.zeros((n_permutations, n_samples))
    
    for i in range(n_permutations):
        sample_order = np.random.permutation(n_samples)
        accumulated_genera = np.zeros(n_samples)
        
        genera_present = set()
        for j, sample in enumerate(sample_order):
            new_genera = set(np.where(genera_matrix[:, sample])[0])
            genera_present.update(new_genera)
            accumulated_genera[j] = len(genera_present)
            
        accumulation_curves[i] = accumulated_genera
    
    mean_curve = np.mean(accumulation_curves, axis=0)
    std_curve = np.std(accumulation_curves, axis=0)
    return mean_curve, std_curve

mean_curve, std_curve = accumulate_genera(genera_matrix)
x = range(1, len(mean_curve) + 1)

plt.figure(figsize=(10, 6))
plt.plot(x, mean_curve, 'b-', label='Mean')
plt.fill_between(x, mean_curve - std_curve, mean_curve + std_curve, alpha=0.2)
plt.xlabel('Number of Samples')
plt.ylabel('Number of Genera')
plt.title('Genera Accumulation Curve in Guma\n(excluding nests S and T)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/genera_accumulation_curve.png', dpi=300, bbox_inches='tight')
plt.show()

# 2. PCoA of genera composition
genera_distances = squareform(pdist(genera_matrix.T, 'jaccard'))
distance_matrix = DistanceMatrix(genera_distances)
pcoa_results = ordination.pcoa(distance_matrix)

# Convert to numpy array for plotting
pcoa_coords = pcoa_results.samples.values

plt.figure(figsize=(10, 8))
unique_nests = np.unique(filtered_nests)
colors = sns.color_palette("husl", n_colors=len(unique_nests))

for nest, color in zip(unique_nests, colors):
    mask = filtered_nests == nest
    mask_indices = np.where(mask)[0]
    
    # Plot points
    plt.scatter(pcoa_coords[mask_indices, 0], 
                pcoa_coords[mask_indices, 1],
                c=[color], 
                label=nest, 
                alpha=0.7)
    
    # Add centroid
    centroid = np.mean(pcoa_coords[mask_indices, :2], axis=0)
    plt.scatter(centroid[0], centroid[1], 
                c='black', 
                marker='x', 
                s=100, 
                linewidth=2)

plt.xlabel(f'PCoA 1 ({pcoa_results.proportion_explained[0]:.1%} explained)')
plt.ylabel(f'PCoA 2 ({pcoa_results.proportion_explained[1]:.1%} explained)')
plt.title('PCoA of Genera Composition in Guma\n(Jaccard dissimilarity)')
plt.legend(title='Nest', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/genera_pcoa.png', dpi=300, bbox_inches='tight')
plt.show()

# 3. Heatmap of genera frequency by nest
genera_freq_by_nest = {}
for nest in unique_nests:
    nest_mask = filtered_nests == nest
    nest_samples = genera_matrix[:, nest_mask]
    genera_freq_by_nest[nest] = np.mean(nest_samples, axis=1) * 100

heatmap_data = pd.DataFrame(genera_freq_by_nest, index=unique_genera)

# Select most frequent genera (top 15)
mean_freq = np.mean(genera_matrix, axis=1)
top_indices = np.argsort(mean_freq)[-15:]
top_genera = [unique_genera[i] for i in top_indices]
heatmap_data_filtered = heatmap_data.loc[top_genera]

plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data_filtered, 
            cmap='YlOrRd', 
            annot=True, 
            fmt='.1f', 
            cbar_kws={'label': 'Detection Frequency (%)'})
plt.title('Top 15 Genera Detection Frequency by Nest in Guma')
plt.ylabel('Genus')
plt.xlabel('Nest')
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/genera_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print(f"Total number of genera: {len(unique_genera)}")
print(f"Average number of genera per sample: {np.mean(np.sum(genera_matrix, axis=0)):.1f}")
print(f"Maximum number of genera in a single sample: {np.max(np.sum(genera_matrix, axis=0))}")
print(f"Minimum number of genera in a single sample: {np.min(np.sum(genera_matrix, axis=0))}")

# Calculate and print genera shared between nests
print("\nNumber of genera shared between nests:")
for i, nest1 in enumerate(unique_nests[:-1]):
    for nest2 in unique_nests[i+1:]:
        nest1_mask = filtered_nests == nest1
        nest2_mask = filtered_nests == nest2
        
        nest1_genera = set(np.where(np.any(genera_matrix[:, nest1_mask], axis=1))[0])
        nest2_genera = set(np.where(np.any(genera_matrix[:, nest2_mask], axis=1))[0])
        
        shared_genera = nest1_genera.intersection(nest2_genera)
        print(f"{nest1} and {nest2}: {len(shared_genera)} genera")