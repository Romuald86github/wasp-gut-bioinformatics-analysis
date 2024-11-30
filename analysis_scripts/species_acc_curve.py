import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Load the presence-absence data
data_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/presence_absence_IK_data.xlsx'
data_df = pd.read_excel(data_path)

# Extract metadata (Nest and Location) and species data
metadata_df = data_df.iloc[0:2, 1:]  # First two rows contain Nest and Location
species_df = data_df.iloc[2:, 1:]    # Species data starts from row 3
species_names = data_df.iloc[2:, 0]   # Species names in first column

# Function to check if a name is a species (has genus and species name)
def is_species_name(name):
    # Split the name and check if it has exactly two parts
    parts = name.split()
    if len(parts) == 2:
        # Check if first part starts with capital letter (genus) and second part is lowercase (species)
        return parts[0][0].isupper() and parts[1].islower()
    return False

# Filter to keep only true species names
species_mask = species_names.apply(is_species_name)
species_df = species_df[species_mask]
species_names = species_names[species_mask]

print(f"Number of species considered: {len(species_names)}")
print("\nFirst few species names:")
print(species_names.head())

# Convert metadata to dictionary format
nests = metadata_df.iloc[0].values
locations = metadata_df.iloc[1].values

# Convert species data to numeric matrix and ensure binary presence/absence
abundance_matrix = species_df.values.astype(float)
abundance_matrix = np.where(abundance_matrix != 0, 1, 0)  # Convert all non-zero values to 1

# Define Hill number function for q=1
def hill_number_q1(abundances):
    """
    Calculate Hill number for q=1 (exponential of Shannon entropy)
    Formula: exp(-sum(p_i * ln(p_i))) where p_i are the relative abundances
    For presence/absence data, p_i = 1/S for all present species, where S is total number of present species
    """
    # Count number of present species
    S = np.sum(abundances > 0)
    if S == 0:
        return 0
        
    # For presence/absence data, all present species have equal relative abundance (1/S)
    p = 1/S
    
    # Calculate Shannon entropy and return its exponential
    # For q=1: H = -sum(p_i * ln(p_i)) = -S * (1/S * ln(1/S)) = ln(S)
    # Therefore, exp(H) = S
    # This means for presence/absence data, Hill number q=1 equals the number of present species
    return S

# Accumulate diversity by location using Hill numbers (q=1)
def accumulate_diversity_by_location(abundance_matrix, locations):
    results = {}
    unique_locations = np.unique(locations)
    
    for location in unique_locations:
        loc_indices = np.where(locations == location)[0]
        loc_abundances = abundance_matrix[:, loc_indices]
        
        # Perform multiple random permutations of the samples
        n_permutations = 100
        temp_accumulation = np.zeros((n_permutations, len(loc_indices)))
        
        for perm in range(n_permutations):
            # Randomly permute the samples
            perm_indices = np.random.permutation(len(loc_indices))
            permuted_data = loc_abundances[:, perm_indices]
            
            # Calculate accumulation
            for i in range(1, permuted_data.shape[1] + 1):
                # Combine samples up to index i
                combined_data = np.any(permuted_data[:, :i] > 0, axis=1).astype(int)
                # Calculate Hill number (q=1)
                diversity = hill_number_q1(combined_data)
                temp_accumulation[perm, i-1] = diversity
        
        # Calculate mean and standard deviation across permutations
        mean_diversity = np.mean(temp_accumulation, axis=0)
        std_diversity = np.std(temp_accumulation, axis=0)
        
        results[location] = {
            'mean': mean_diversity,
            'std': std_diversity
        }
    
    return results

# Calculate diversity accumulation
diversity_accumulation = accumulate_diversity_by_location(abundance_matrix, locations)

# Plot the diversity accumulation curve
def plot_diversity_accumulation(diversity_data, save_dir):
    plt.figure(figsize=(10, 7))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(diversity_data)))
    
    for (location, values), color in zip(diversity_data.items(), colors):
        x = np.arange(1, len(values['mean']) + 1)
        
        if len(x) > 3:
            # Smoothing using spline interpolation
            x_smooth = np.linspace(x.min(), x.max(), 300)
            
            # Interpolate both mean and standard deviation
            spline_mean = make_interp_spline(x, values['mean'], k=min(3, len(x)-1))
            spline_std = make_interp_spline(x, values['std'], k=min(3, len(x)-1))
            
            y_smooth = spline_mean(x_smooth)
            std_smooth = spline_std(x_smooth)
            
            # Calculate confidence intervals
            lower_bound = y_smooth - std_smooth
            upper_bound = y_smooth + std_smooth
            
            # Plot smoothed curve and confidence interval
            plt.plot(x_smooth, y_smooth, label=f"{location} (n={len(x)})", color=color)
            plt.fill_between(x_smooth, lower_bound, upper_bound, alpha=0.2, color=color)
        else:
            # For locations with few points, plot direct lines
            plt.plot(x, values['mean'], label=f"{location} (n={len(x)})", marker='o', color=color)
            plt.fill_between(x, values['mean'] - values['std'], 
                           values['mean'] + values['std'], alpha=0.2, color=color)
    
    plt.title(f"Species Diversity Accumulation Curve by Location\n(Hill number q=1, Total Species: {len(species_names)})")
    plt.xlabel("Number of Samples")
    plt.ylabel("Effective Number of Species (Hill q=1)")
    plt.legend(title="Location", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(f'{save_dir}/diversity_accumulation_curve_hill_q1.png', bbox_inches='tight', dpi=300)
    plt.show()

# Define the save directory path
save_dir = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics'

# Plot and save the curve
plot_diversity_accumulation(diversity_accumulation, save_dir)