import pandas as pd
import sys

# Define default paths
merged_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results.tsv'
sample_origin_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/sample_origin_mapping_cleaned.tsv'
updated_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'

# Override paths if provided as command line arguments
if len(sys.argv) == 4:
   merged_file_path = sys.argv[1]
   sample_origin_file_path = sys.argv[2]
   updated_file_path = sys.argv[3]

# Load and process data
merged_df = pd.read_csv(merged_file_path, sep='\t', index_col=0)
sample_origin_df = pd.read_csv(sample_origin_file_path, sep='\t')

# Create mapping and rename columns
sample_to_origin = dict(zip(sample_origin_df['sample'], sample_origin_df['origin']))
merged_df.rename(columns=sample_to_origin, inplace=True)

# Save results
merged_df.to_csv(updated_file_path, sep='\t', index=False)
print(f"merged table with origins saved as: {updated_file_path}")