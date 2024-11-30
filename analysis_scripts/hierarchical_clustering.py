import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import networkx as nx
from collections import Counter

# Load and process data (same as before)
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

# Extract genera and create matrices (same as before)
def get_genus(name):
    return name.split()[0]

genera_names = original_names.apply(get_genus)
species_names = original_names

# Create genera presence/absence matrix
genera_to_rows = {}
for idx, genus in enumerate(genera_names):
    if genus not in genera_to_rows:
        genera_to_rows[genus] = []
    genera_to_rows[genus].append(idx)

abundance_matrix = filtered_species_df.values.astype(float)
abundance_matrix = np.where(abundance_matrix != 0, 1, 0)

unique_genera = sorted(genera_to_rows.keys())
genera_matrix = np.zeros((len(unique_genera), abundance_matrix.shape[1]))
for i, genus in enumerate(unique_genera):
    row_indices = genera_to_rows[genus]
    genera_matrix[i, :] = np.any(abundance_matrix[row_indices, :], axis=0)

# 1. Species and Genera Richness Boxplot by Nest
species_richness = np.sum(abundance_matrix, axis=0)
genera_richness = np.sum(genera_matrix, axis=0)

richness_df = pd.DataFrame({
    'Nest': np.repeat(filtered_nests, 2),
    'Richness': np.concatenate([species_richness, genera_richness]),
    'Level': ['Species'] * len(filtered_nests) + ['Genus'] * len(filtered_nests)
})

plt.figure(figsize=(12, 6))
sns.boxplot(data=richness_df, x='Nest', y='Richness', hue='Level')
plt.title('Species and Genera Richness by Nest')
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/richness_comparison.png', dpi=300)
plt.show()

# Perform Kruskal-Wallis test for richness differences between nests
print("\nKruskal-Wallis Test Results:")
print("Species Richness between nests:")
h_stat, p_val = kruskal(*[species_richness[filtered_nests == nest] for nest in np.unique(filtered_nests)])
print(f"H-statistic: {h_stat:.3f}, p-value: {p_val:.3f}")

print("\nGenera Richness between nests:")
h_stat, p_val = kruskal(*[genera_richness[filtered_nests == nest] for nest in np.unique(filtered_nests)])
print(f"H-statistic: {h_stat:.3f}, p-value: {p_val:.3f}")

# 2. Hierarchical Clustering of Samples
# Calculate distance matrix using Jaccard distance
dist_matrix = pdist(genera_matrix.T, 'jaccard')
linkage_matrix = linkage(dist_matrix, method='ward')

plt.figure(figsize=(12, 6))
dendrogram(linkage_matrix, labels=filtered_nests, leaf_rotation=90)
plt.title('Hierarchical Clustering of Samples Based on Genera Composition')
plt.xlabel('Sample (Nest)')
plt.ylabel('Distance')
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/hierarchical_clustering.png', dpi=300)
plt.show()

# 3. Species-Genera Ratio Analysis
species_per_genus = pd.Series(genera_names).value_counts()
plt.figure(figsize=(10, 6))
species_per_genus.head(15).plot(kind='bar')
plt.title('Number of Species per Genus (Top 15)')
plt.xlabel('Genus')
plt.ylabel('Number of Species')
plt.xticks(rotation=45, ha='right')
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/species_per_genus.png', dpi=300)
plt.show()

# 4. Co-occurrence Network (for top genera)
# Calculate co-occurrence matrix
top_genera_idx = np.argsort(np.sum(genera_matrix, axis=1))[-10:]  # top 10 genera
co_occurrence = np.zeros((len(top_genera_idx), len(top_genera_idx)))

for i, idx1 in enumerate(top_genera_idx):
    for j, idx2 in enumerate(top_genera_idx):
        if i < j:
            co_occur = np.sum(np.logical_and(genera_matrix[idx1, :], genera_matrix[idx2, :]))
            co_occurrence[i, j] = co_occur
            co_occurrence[j, i] = co_occur

plt.figure(figsize=(10, 10))
G = nx.Graph()
genera_names = [unique_genera[i] for i in top_genera_idx]

# Add nodes
for genus in genera_names:
    G.add_node(genus)

# Add edges
for i in range(len(genera_names)):
    for j in range(i+1, len(genera_names)):
        if co_occurrence[i, j] > 0:
            G.add_edge(genera_names[i], genera_names[j], weight=co_occurrence[i, j])

# Draw network
pos = nx.spring_layout(G)
edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]
nx.draw(G, pos, 
        with_labels=True,
        node_color='lightblue',
        node_size=1000,
        font_size=8,
        width=[w/2 for w in weights],
        edge_color='gray')
plt.title('Co-occurrence Network of Top 10 Genera')
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/cooccurrence_network.png', dpi=300)
plt.show()

# 5. Additional Statistics
print("\nAdditional Statistics:")
print(f"\nSpecies per Genus Statistics:")
print(f"Mean species per genus: {np.mean(species_per_genus):.2f}")
print(f"Median species per genus: {np.median(species_per_genus):.2f}")
print(f"Max species in a genus: {np.max(species_per_genus)}")
print(f"Genera with single species: {np.sum(species_per_genus == 1)}")

# Calculate sample coverage
sample_coverage = np.sum(genera_matrix, axis=0)
print(f"\nSample Coverage Statistics:")
print(f"Mean genera per sample: {np.mean(sample_coverage):.2f}")
print(f"Standard deviation: {np.std(sample_coverage):.2f}")
print(f"Coefficient of variation: {np.std(sample_coverage)/np.mean(sample_coverage):.2f}")

# Unique combinations analysis
print("\nUnique Combinations Analysis:")
unique_combinations = {}
for nest in np.unique(filtered_nests):
    nest_mask = filtered_nests == nest
    nest_genera = genera_matrix[:, nest_mask]
    unique_patterns = set(tuple(col) for col in nest_genera.T)
    unique_combinations[nest] = len(unique_patterns)
    print(f"Nest {nest}: {len(unique_patterns)} unique genera combinations")