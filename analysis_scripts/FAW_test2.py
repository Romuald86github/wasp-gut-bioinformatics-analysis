import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

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

# Calculate observed detection rates per nest
nest_detection_rates = {}
nest_sample_sizes = {}
for nest in unique_nests:
    nest_mask = filtered_nests == nest
    nest_detection_rates[nest] = np.mean(spodoptera_presence[nest_mask])
    nest_sample_sizes[nest] = np.sum(nest_mask)

# Calculate overall detection rate (this will be our null expectation for each nest)
overall_detection_rate = np.mean(spodoptera_presence)

# Randomization test
n_permutations = 10000

# For each nest, calculate how likely the observed detection rate is under equal detection rates
nest_random_rates = {nest: np.zeros(n_permutations) for nest in unique_nests}
for nest in unique_nests:
    n_samples = nest_sample_sizes[nest]
    for i in range(n_permutations):
        # Generate random samples with overall detection rate
        random_samples = np.random.random(n_samples) < overall_detection_rate
        nest_random_rates[nest][i] = np.mean(random_samples)

# Calculate two-tailed p-values for each nest
nest_p_values = {}
for nest in unique_nests:
    observed = nest_detection_rates[nest]
    random_rates = nest_random_rates[nest]
    p_value = 2 * min(
        np.mean(random_rates >= observed),
        np.mean(random_rates <= observed)
    )
    nest_p_values[nest] = p_value

# Visualization
plt.figure(figsize=(12, 6))

# Plot nest-specific detection rates vs null expectation
nest_names = list(nest_detection_rates.keys())
observed_rates = [nest_detection_rates[nest] for nest in nest_names]
random_means = [np.mean(nest_random_rates[nest]) for nest in nest_names]
random_std = [np.std(nest_random_rates[nest]) for nest in nest_names]

x = range(len(nest_names))
plt.errorbar(x, [overall_detection_rate] * len(nest_names), 
            yerr=[2*np.std(nest_random_rates[nest]) for nest in nest_names], 
            fmt='o', label='Null expectation (mean ± 2SD)',
            color='blue', alpha=0.5)
plt.scatter(x, observed_rates, color='red', label='Observed', zorder=5, s=100)
plt.axhline(y=overall_detection_rate, color='blue', linestyle='--', 
           label='Overall detection rate')
plt.xticks(x, nest_names, rotation=45)
plt.xlabel('Nest')
plt.ylabel('Detection Rate')
plt.title('Observed vs Expected Detection Rates by Nest\n(Null: equal detection rates)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/spodoptera_equal_rates_test.png', 
            dpi=300, bbox_inches='tight')
plt.show()

# Print results
print("Test Results for Equal Detection Rates of S. frugiperda")
print("-" * 50)
print(f"\nOverall detection rate (null expectation): {overall_detection_rate:.3f}")

print("\nNest-specific Results:")
print("-" * 30)
for nest in nest_names:
    obs_rate = nest_detection_rates[nest]
    n_samples = nest_sample_sizes[nest]
    n_detections = int(obs_rate * n_samples)
    exp_detections = overall_detection_rate * n_samples
    p_val = nest_p_values[nest]
    
    print(f"\nNest {nest}:")
    print(f"Samples: {n_samples}")
    print(f"Observed detections: {n_detections} ({obs_rate:.3f})")
    print(f"Expected detections: {exp_detections:.1f} ({overall_detection_rate:.3f})")
    print(f"P-value (two-tailed): {p_val:.4f}")
    if p_val < 0.05:
        if obs_rate > overall_detection_rate:
            print("* Significantly higher than expected *")
        else:
            print("* Significantly lower than expected *")

# Chi-square test for homogeneity across nests
observed_counts = []
expected_counts = []
for nest in unique_nests:
    nest_mask = filtered_nests == nest
    n_samples = np.sum(nest_mask)
    observed_counts.append(np.sum(spodoptera_presence[nest_mask]))
    expected_counts.append(overall_detection_rate * n_samples)

chi2, chi2_p = stats.chisquare(observed_counts, expected_counts)

print("\nChi-square Test for Homogeneity:")
print("-" * 30)
print(f"Chi-square statistic: {chi2:.3f}")
print(f"P-value: {chi2_p:.4f}")

# Calculate effect sizes (deviation from expected rate)
print("\nEffect Sizes (deviation from expected rate):")
print("-" * 30)
for nest in unique_nests:
    obs_rate = nest_detection_rates[nest]
    effect_size = (obs_rate - overall_detection_rate) / overall_detection_rate * 100
    print(f"Nest {nest}: {effect_size:+.1f}% deviation from expected")