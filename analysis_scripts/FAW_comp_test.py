import pandas as pd
import numpy as np
from scipy import stats
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

# Find Spodoptera frugiperda
spodoptera_mask = original_names.str.contains('Spodoptera frugiperda', case=False)
spodoptera_index = np.where(spodoptera_mask)[0]
if len(spodoptera_index) == 0:
    raise ValueError("Spodoptera frugiperda not found in species list")

spodoptera_presence = abundance_matrix[spodoptera_index[0], :]

# Calculate observed statistics
unique_nests = np.unique(filtered_nests)
n_nests = len(unique_nests)

# Calculate observed detection counts and rates per nest
nest_detection_counts = {}
nest_sample_sizes = {}
nest_detection_rates = {}
for nest in unique_nests:
    nest_mask = filtered_nests == nest
    nest_detection_counts[nest] = np.sum(spodoptera_presence[nest_mask])
    nest_sample_sizes[nest] = np.sum(nest_mask)
    nest_detection_rates[nest] = nest_detection_counts[nest] / nest_sample_sizes[nest]

# Calculate overall detection rate
total_detections = sum(nest_detection_counts.values())
total_samples = sum(nest_sample_sizes.values())
overall_detection_rate = total_detections / total_samples

# Perform randomization test
n_permutations = 10000
observed_stat = np.sum([(rate - overall_detection_rate)**2 for rate in nest_detection_rates.values()])

# Calculate random stats under null hypothesis of equal detection rates
random_stats = np.zeros(n_permutations)
for i in range(n_permutations):
    random_presence = np.random.binomial(1, overall_detection_rate, total_samples)
    random_rates = {}
    start_idx = 0
    for nest in unique_nests:
        n_samples = nest_sample_sizes[nest]
        end_idx = start_idx + n_samples
        random_rates[nest] = np.mean(random_presence[start_idx:end_idx])
        start_idx = end_idx
    random_stats[i] = np.sum([(rate - overall_detection_rate)**2 for rate in random_rates.values()])

# Calculate p-value
p_value = np.mean(random_stats >= observed_stat)

# Visualization
plt.figure(figsize=(12, 6))

# Plot 1: Detection rates by nest
plt.subplot(1, 2, 1)
x = range(len(unique_nests))
plt.scatter(x, [nest_detection_rates[nest] for nest in unique_nests], 
           color='red', s=100, zorder=5)
plt.axhline(y=overall_detection_rate, color='blue', linestyle='--', 
           label='Expected rate')
plt.fill_between([-0.5, len(unique_nests)-0.5], 
                 [overall_detection_rate - 2*np.std(list(nest_detection_rates.values()))]*2,
                 [overall_detection_rate + 2*np.std(list(nest_detection_rates.values()))]*2,
                 color='blue', alpha=0.2, label='±2 SD range')
plt.xticks(x, unique_nests, rotation=45)
plt.xlabel('Nest')
plt.ylabel('Detection Rate')
plt.title('S. frugiperda Detection Rates by Nest')
plt.legend()

# Plot 2: Null distribution
plt.subplot(1, 2, 2)
plt.hist(random_stats, bins=50, alpha=0.7, density=True)
plt.axvline(observed_stat, color='red', linestyle='--',
            label=f'Observed (p={p_value:.4f})')
plt.xlabel('Test Statistic')
plt.ylabel('Density')
plt.title('Null Distribution')
plt.legend()

plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/spodoptera_detection_test.png', 
            dpi=300, bbox_inches='tight')
plt.show()

# Print results
print("Test for Equal Detection Rates of S. frugiperda Across Nests")
print("-" * 60)
print(f"\nOverall detection rate: {overall_detection_rate:.3f}")
print(f"Total detections: {total_detections} out of {total_samples} samples")
print(f"\nOverall test p-value: {p_value:.4f}")
if p_value < 0.05:
    print("* Detection rates are significantly different across nests *")
else:
    print("* No significant evidence of different detection rates across nests *")

print("\nDetection rates by nest:")
print("-" * 30)
for nest in unique_nests:
    print(f"Nest {nest}: {nest_detection_counts[nest]} out of {nest_sample_sizes[nest]} " +
          f"({nest_detection_rates[nest]:.3f})")