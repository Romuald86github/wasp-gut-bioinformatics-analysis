import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

# Create presence/absence matrix
abundance_matrix = filtered_species_df.values.astype(float)
abundance_matrix = np.where(abundance_matrix != 0, 1, 0)

# List of target species
target_species = [
    'Spodoptera frugiperda',
    'Thaumatotibia leucotreta',
    'Spodoptera littoralis',
    'Maruca vitrata',
    'Helicoverpa assulta'
]

# Common names for labels
common_names = {
    'Spodoptera frugiperda': 'Fall armyworm',
    'Thaumatotibia leucotreta': 'False codling moth',
    'Spodoptera littoralis': 'African cotton leafworm',
    'Maruca vitrata': 'Legume pod borer',
    'Helicoverpa assulta': 'Oriental tobacco budworm'
}

# Get detection data for each species
species_data = {}
for species in target_species:
    species_mask = original_names.str.contains(species, case=False)
    if any(species_mask):
        species_index = np.where(species_mask)[0][0]
        species_presence = abundance_matrix[species_index, :]
        
        # Overall detection rate
        overall_rate = np.mean(species_presence)
        
        # Detection rates by nest
        nest_rates = {}
        for nest in np.unique(filtered_nests):
            nest_mask = filtered_nests == nest
            nest_rates[nest] = np.mean(species_presence[nest_mask])
        
        species_data[species] = {
            'overall_rate': overall_rate,
            'nest_rates': nest_rates,
            'presence_data': species_presence
        }

# Create plots
plt.figure(figsize=(15, 10))

# Plot 1: Overall detection rates
plt.subplot(2, 1, 1)
species_names = list(species_data.keys())
overall_rates = [species_data[sp]['overall_rate'] * 100 for sp in species_names]
labels = [f"{name}\n({common_names[name]})" for name in species_names]

bars = plt.bar(range(len(species_names)), overall_rates)
plt.xticks(range(len(species_names)), labels, rotation=45, ha='right')
plt.ylabel('Detection Rate (%)')
plt.title('Overall Detection Rates of Key Pest Species in Guma')

# Add percentage labels on top of bars
for i, bar in enumerate(bars):
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{overall_rates[i]:.1f}%',
             ha='center', va='bottom')

# Plot 2: Detection rates by nest
plt.subplot(2, 1, 2)
data_for_heatmap = np.array([[species_data[sp]['nest_rates'][nest] * 100 
                             for nest in np.unique(filtered_nests)] 
                            for sp in species_names])

sns.heatmap(data_for_heatmap,
            annot=True,
            fmt='.1f',
            cmap='YlOrRd',
            xticklabels=np.unique(filtered_nests),
            yticklabels=labels,
            cbar_kws={'label': 'Detection Rate (%)'})
plt.xlabel('Nest')
plt.title('Detection Rates by Nest')

plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/pest_detection_comparison.png', 
            dpi=300, bbox_inches='tight')
plt.show()

# Print detailed statistics
print("Detection Statistics for Key Pest Species in Guma")
print("-" * 60)

for species in species_names:
    data = species_data[species]
    total_detections = np.sum(data['presence_data'])
    total_samples = len(data['presence_data'])
    
    print(f"\n{species} ({common_names[species]}):")
    print(f"Overall detection rate: {data['overall_rate']*100:.1f}% ({total_detections} out of {total_samples} samples)")
    print("\nDetection rates by nest:")
    for nest, rate in data['nest_rates'].items():
        nest_mask = filtered_nests == nest
        nest_detections = np.sum(data['presence_data'][nest_mask])
        nest_samples = np.sum(nest_mask)
        print(f"  Nest {nest}: {rate*100:.1f}% ({nest_detections} out of {nest_samples} samples)")

# Perform chi-square test for homogeneity
print("\nChi-square test for homogeneity of detection rates:")
print("-" * 60)
observed = np.array([np.sum(species_data[sp]['presence_data']) for sp in species_names])
expected = np.mean(observed) * np.ones_like(observed)
chi2, p_value = stats.chisquare(observed, expected)
print(f"Chi-square statistic: {chi2:.2f}")
print(f"P-value: {p_value:.4f}")
if p_value < 0.05:
    print("* Detection rates are significantly different among species *")
else:
    print("* No significant evidence of different detection rates among species *")