import os
import pandas as pd
import re
import sys

# Define default paths
cutadapt_folder_path = '/Users/romualdchristialtcheutchoua/Downloads/IK_data/cutadapt'
output_file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/sample_origin_mapping_cleaned.tsv'

# Override paths if provided as command line arguments
if len(sys.argv) == 3:
   cutadapt_folder_path = sys.argv[1]
   output_file_path = sys.argv[2]

# Rest of the script remains the same
pattern = re.compile(r'([A-Z]\d+)_S(\d+)_L001_R1_001.fastq.gz')
samples = []
origins = []

for filename in os.listdir(cutadapt_folder_path):
   match = pattern.search(filename)
   if match:
       origin = match.group(1)
       sample = f"S{match.group(2)}"
       samples.append(sample)
       origins.append(origin)

sample_origin_df = pd.DataFrame({
   'sample': samples,
   'origin': origins
})

sample_origin_df.sort_values('sample', inplace=True)

print("Duplicate samples and their origins:")
duplicates = sample_origin_df[sample_origin_df.duplicated('origin', keep=False)]
for index, row in duplicates.iterrows():
   print(f"Sample: {row['sample']}, Origin: {row['origin']}")

sample_origin_df.drop_duplicates('origin', keep='first', inplace=True)

print("\nCleaned DataFrame:")
print(sample_origin_df)

sample_origin_df.to_csv(output_file_path, sep='\t', index=False)
print(f"\nCleaned sample-origin mapping saved as: {output_file_path}")