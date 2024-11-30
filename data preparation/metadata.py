import pandas as pd

# Load the origin data (which contains sample origins, nests, and larvae)
origin_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/wasp_gut_results_with_origins.tsv'
origin_df = pd.read_csv(origin_path, sep='\t')

# (1) Remove the ASV column (first column) and focus on the remaining columns
origin_df = origin_df.drop(columns=['ASV'])

# (2) Extract the larvae samples (columns excluding the last 8 taxonomy columns)
taxonomy_cols = ['Order', 'Family', 'Genus', 'Species']
sample_cols = [col for col in origin_df.columns if col not in taxonomy_cols]

# Create a metadata DataFrame from the sample columns
metadata = pd.DataFrame()

# Extract nest IDs (first character of each sample column)
metadata['Nest'] = [sample[0] for sample in sample_cols]

# Extract larvae IDs (full column names)
metadata['Larvae'] = sample_cols

# Define locations based on nest IDs
def assign_location(nest_id):
    if nest_id in 'ABCDEFGHIJKLMNOP':  # Nests A-P sampled in Guma
        return 'Guma'
    elif nest_id == 'Q':  # Nest Q sampled in Makurdi
        return 'Makurdi'
    elif nest_id in 'RST':  # Nests R-T sampled in Guma
        return 'Guma'
    elif nest_id == 'U':  # Nest U sampled in Nsukka
        return 'Nsukka'
    elif nest_id == 'W':  # Nest W sampled in Makurdi
        return 'Makurdi'
    elif nest_id == 'X':  # Nest X sampled in Guma
        return 'Guma'
    else:
        return 'Unknown'

metadata['Location'] = metadata['Nest'].apply(assign_location)

# Set 'Larvae' as the index
metadata.set_index('Larvae', inplace=True)

# Save the metadata to the specified path as 'metadata.tsv'
metadata_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/metadata.tsv'
metadata.to_csv(metadata_path, sep='\t')

# (3) Now that we have the metadata (Nest, Larvae, Location), let's print it
print(metadata.head())
