import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.spatial.distance import pdist, squareform

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

# Extract genera
def get_genus(name):
    return name.split()[0]

genera_names = original_names.apply(get_genus)
species_names = original_names

# Create matrices
abundance_matrix = filtered_species_df.values.astype(float)
abundance_matrix = np.where(abundance_matrix != 0, 1, 0)

# Create genera matrix
genera_to_rows = {}
for idx, genus in enumerate(genera_names):
    if genus not in genera_to_rows:
        genera_to_rows[genus] = []
    genera_to_rows[genus].append(idx)

unique_genera = sorted(genera_to_rows.keys())
genera_matrix = np.zeros((len(unique_genera), abundance_matrix.shape[1]))
for i, genus in enumerate(unique_genera):
    row_indices = genera_to_rows[genus]
    genera_matrix[i, :] = np.any(abundance_matrix[row_indices, :], axis=0)

# 1. Network of Nests based on Shared Species
def create_nest_network(abundance_data, nests, threshold=0.3):
    unique_nests = np.unique(nests)
    G = nx.Graph()
    
    # Add nodes (nests)
    for nest in unique_nests:
        G.add_node(nest, node_type='nest')
        
    # Calculate shared species/genera between nests
    for i, nest1 in enumerate(unique_nests):
        for j, nest2 in enumerate(unique_nests):
            if i < j:
                mask1 = nests == nest1
                mask2 = nests == nest2
                
                # Get species/genera present in each nest
                nest1_present = np.any(abundance_data[:, mask1], axis=1)
                nest2_present = np.any(abundance_data[:, mask2], axis=1)
                
                # Calculate Jaccard similarity
                intersection = np.sum(nest1_present & nest2_present)
                union = np.sum(nest1_present | nest2_present)
                jaccard = intersection / union if union > 0 else 0
                
                if jaccard >= threshold:
                    G.add_edge(nest1, nest2, weight=jaccard, shared=intersection)
    
    return G

# Create and plot species network
plt.figure(figsize=(12, 8))
G_species = create_nest_network(abundance_matrix, filtered_nests, threshold=0.3)

pos = nx.spring_layout(G_species, k=1, iterations=50)
nx.draw_networkx_nodes(G_species, pos, node_color='lightblue', 
                      node_size=1000, alpha=0.6)
nx.draw_networkx_labels(G_species, pos)

# Draw edges with width proportional to Jaccard similarity
edges = G_species.edges(data=True)
weights = [d['weight']*3 for (u, v, d) in edges]
nx.draw_networkx_edges(G_species, pos, width=weights, alpha=0.5)

# Add edge labels showing number of shared species
edge_labels = {(u, v): f'{d["shared"]} sp.' for (u, v, d) in edges}
nx.draw_networkx_edge_labels(G_species, pos, edge_labels, font_size=8)

plt.title('Network of Nests Connected by Shared Species\n(Edge width proportional to Jaccard similarity)')
plt.axis('off')
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/nest_species_network.png', dpi=300, bbox_inches='tight')
plt.show()

# Create and plot genera network
plt.figure(figsize=(12, 8))
G_genera = create_nest_network(genera_matrix, filtered_nests, threshold=0.3)

pos = nx.spring_layout(G_genera, k=1, iterations=50)
nx.draw_networkx_nodes(G_genera, pos, node_color='lightgreen', 
                      node_size=1000, alpha=0.6)
nx.draw_networkx_labels(G_genera, pos)

# Draw edges with width proportional to Jaccard similarity
edges = G_genera.edges(data=True)
weights = [d['weight']*3 for (u, v, d) in edges]
nx.draw_networkx_edges(G_genera, pos, width=weights, alpha=0.5)

# Add edge labels showing number of shared genera
edge_labels = {(u, v): f'{d["shared"]} gen.' for (u, v, d) in edges}
nx.draw_networkx_edge_labels(G_genera, pos, edge_labels, font_size=8)

plt.title('Network of Nests Connected by Shared Genera\n(Edge width proportional to Jaccard similarity)')
plt.axis('off')
plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/nest_genera_network.png', dpi=300, bbox_inches='tight')
plt.show()

# 2. Print detailed information about shared species/genera
print("\nDetailed Analysis of Shared Taxa between Nests:")
unique_nests = np.unique(filtered_nests)
for i, nest1 in enumerate(unique_nests):
    for j, nest2 in enumerate(unique_nests):
        if i < j:
            mask1 = filtered_nests == nest1
            mask2 = filtered_nests == nest2
            
            # Species analysis
            species1 = set(np.where(np.any(abundance_matrix[:, mask1], axis=1))[0])
            species2 = set(np.where(np.any(abundance_matrix[:, mask2], axis=1))[0])
            shared_species = species1.intersection(species2)
            
            # Genera analysis
            genera1 = set(np.where(np.any(genera_matrix[:, mask1], axis=1))[0])
            genera2 = set(np.where(np.any(genera_matrix[:, mask2], axis=1))[0])
            shared_genera = genera1.intersection(genera2)
            
            print(f"\nNests {nest1} and {nest2}:")
            print(f"Shared species: {len(shared_species)}")
            print(f"Shared genera: {len(shared_genera)}")
            
            # List some shared species as examples
            if shared_species:
                print("Example shared species:")
                for idx in list(shared_species)[:5]:  # Show up to 5 examples
                    print(f"  - {species_names.iloc[idx]}")
            
            # List shared genera
            if shared_genera:
                print("Shared genera:")
                for idx in shared_genera:
                    print(f"  - {unique_genera[idx]}")

# 3. Calculate and print network metrics
print("\nNetwork Metrics:")

print("\nSpecies Network:")
print(f"Average clustering coefficient: {nx.average_clustering(G_species):.3f}")
print(f"Network density: {nx.density(G_species):.3f}")
if nx.is_connected(G_species):
    print(f"Average shortest path length: {nx.average_shortest_path_length(G_species):.3f}")
print(f"Number of connected components: {nx.number_connected_components(G_species)}")

print("\nGenera Network:")
print(f"Average clustering coefficient: {nx.average_clustering(G_genera):.3f}")
print(f"Network density: {nx.density(G_genera):.3f}")
if nx.is_connected(G_genera):
    print(f"Average shortest path length: {nx.average_shortest_path_length(G_genera):.3f}")
print(f"Number of connected components: {nx.number_connected_components(G_genera)}")