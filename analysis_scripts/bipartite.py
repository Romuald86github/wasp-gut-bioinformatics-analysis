import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

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

def create_bipartite_network(data_matrix, nests, taxa_names, title, filename, top_n=15):
    # Create graph
    G = nx.Graph()
    
    # Get top taxa by frequency
    taxa_freq = np.sum(data_matrix, axis=1)
    top_taxa_idx = np.argsort(taxa_freq)[-top_n:]
    
    # Add nest nodes
    nest_positions = {}
    unique_nests = np.unique(nests)
    n_nests = len(unique_nests)
    for i, nest in enumerate(unique_nests):
        G.add_node(nest, bipartite=0)  # bipartite=0 for nests
        nest_positions[nest] = (i * (n_nests/2), 1.0)
    
    # Add taxa nodes and edges
    taxa_positions = {}
    for i, idx in enumerate(top_taxa_idx):
        taxon = taxa_names.iloc[idx] if isinstance(taxa_names, pd.Series) else taxa_names[idx]
        G.add_node(taxon, bipartite=1)  # bipartite=1 for taxa
        taxa_positions[taxon] = (i * (top_n/3), 0.0)
        
        # Add edges
        for j, nest in enumerate(unique_nests):
            nest_mask = nests == nest
            if np.any(data_matrix[idx, nest_mask]):
                G.add_edge(nest, taxon)
    
    # Combine positions
    pos = {**nest_positions, **taxa_positions}
    
    # Create plot
    plt.figure(figsize=(15, 8))
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, 
                          nodelist=[n for n, d in G.nodes(data=True) if d['bipartite']==0],
                          node_color='lightblue', node_size=1000, alpha=0.7)
    nx.draw_networkx_nodes(G, pos,
                          nodelist=[n for n, d in G.nodes(data=True) if d['bipartite']==1],
                          node_color='lightgreen', node_size=1000, alpha=0.7)
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos)
    
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/{filename}.png', 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    return G

# Create and plot species network
G_species = create_bipartite_network(
    abundance_matrix, 
    filtered_nests,
    species_names,
    'Bipartite Network: Nests and Top 15 Most Common Species',
    'nest_species_bipartite'
)

# Create and plot genera network
G_genera = create_bipartite_network(
    genera_matrix,
    filtered_nests,
    unique_genera,
    'Bipartite Network: Nests and Top 15 Most Common Genera',
    'nest_genera_bipartite'
)

# Print network statistics
def print_network_stats(G, network_type):
    print(f"\n{network_type} Network Statistics:")
    
    # Get the two node sets
    top_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==0}  # nests
    bottom_nodes = {n for n, d in G.nodes(data=True) if d['bipartite']==1}  # taxa
    
    print(f"Number of nests: {len(top_nodes)}")
    print(f"Number of taxa: {len(bottom_nodes)}")
    print(f"Total number of connections: {G.number_of_edges()}")
    
    # Calculate average degree for each set
    top_degrees = [d for n, d in G.degree() if n in top_nodes]
    bottom_degrees = [d for n, d in G.degree() if n in bottom_nodes]
    
    print(f"Average connections per nest: {np.mean(top_degrees):.2f}")
    print(f"Average connections per taxon: {np.mean(bottom_degrees):.2f}")
    
    # Print most connected nodes
    print("\nMost connected nests:")
    nest_degrees = [(n, d) for n, d in G.degree() if n in top_nodes]
    for nest, degree in sorted(nest_degrees, key=lambda x: x[1], reverse=True):
        print(f"{nest}: {degree} connections")
    
    print("\nMost connected taxa:")
    taxa_degrees = [(n, d) for n, d in G.degree() if n in bottom_nodes]
    for taxon, degree in sorted(taxa_degrees, key=lambda x: x[1], reverse=True):
        print(f"{taxon}: {degree} connections")

# Print statistics for both networks
print_network_stats(G_species, "Species")
print_network_stats(G_genera, "Genera")