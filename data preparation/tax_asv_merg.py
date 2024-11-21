import pandas as pd
import sys

# Define default file paths
asv_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/table_ASVs.tsv'
taxonomy_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/Final taxonomy .csv'
output_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results.tsv'

# Override paths if provided as command line arguments
if len(sys.argv) == 4:
    asv_file_path = sys.argv[1]
    taxonomy_file_path = sys.argv[2]
    output_file_path = sys.argv[3]

# Load the ASV count table
asv_df = pd.read_csv(asv_file_path, sep='\t')

# Load the taxonomy table, skipping the first row (original header)
taxonomy_df = pd.read_csv(taxonomy_file_path)

# Filter asv_df to include only ASVs present in the taxonomy_df
asv_df_filtered = asv_df[asv_df['ASV'].isin(taxonomy_df['ASV'])]

# Merge the filtered ASV count table with the taxonomy table on the 'ASV' column
merged_df = pd.merge(asv_df_filtered, taxonomy_df, on='ASV', how='left')
merged_df = merged_df.drop(columns=['Similarity percentage'])

# Save the resulting dataframe to a new TSV file
merged_df.to_csv(output_file_path, sep='\t', index=True)
print(f"Merged table saved as: {output_file_path}")