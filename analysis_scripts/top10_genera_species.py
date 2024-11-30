import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
import warnings
warnings.filterwarnings('ignore')

# Load the data
file_path = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics/presence_absence_IK_data.xlsx'
data = pd.read_excel(file_path)

# Extract metadata from first two rows
nests = data.iloc[0, 1:]  # First row contains nest information
locations = data.iloc[1, 1:]  # Second row contains location information

# Clean and prepare the data
cleaned_data = data.drop([0, 1]).reset_index(drop=True)
cleaned_data.columns.values[0] = 'Taxa'

# Filter columns for Guma location excluding nests S and T
guma_columns = [col for col in cleaned_data.columns[1:] 
                if locations[col] == "Guma" and nests[col] not in ["S", "T"]]
guma_data = cleaned_data[['Taxa'] + guma_columns]

# Convert to presence/absence (1/0)
guma_data[guma_columns] = guma_data[guma_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
guma_data[guma_columns] = guma_data[guma_columns].applymap(lambda x: 1 if x > 0 else 0)

# Count number of larvae where each taxon was detected
guma_data['Larvae_Detections'] = guma_data[guma_columns].sum(axis=1)

# Extract genus and identify species
guma_data['Genus'] = guma_data['Taxa'].apply(lambda x: x.split()[0])
guma_data['Is_Species'] = guma_data['Taxa'].apply(lambda x: len(x.split()) > 1)

# Get top 10 species by number of larvae detections
species_data = (guma_data[guma_data['Is_Species']]
               .sort_values(by='Larvae_Detections', ascending=False)
               .head(10))

# Calculate genus detections (number of larvae where any species of the genus was detected)
genus_detections = {}
for genus in guma_data['Genus'].unique():
    genus_mask = guma_data['Genus'] == genus
    genus_samples = guma_data[genus_mask][guma_columns].max()  # Use max to get presence in any species
    genus_detections[genus] = genus_samples.sum()

# Convert to series and get top 10
genera_data = (pd.Series(genus_detections)
              .sort_values(ascending=False)
              .head(10))

# Create directory path
save_dir = '/Users/romualdchristialtcheutchoua/Desktop/IK_bioinformatics'

# Create simple SVG icons for different genera (you might want to customize these)
def create_moth_icon(style='default'):
    icons = {
        'default': '''
        <svg viewBox="0 0 50 50">
            <path fill="black" d="
                M25,10 L35,25 L25,40 L15,25 Z
                M25,20 L30,25 L25,30 L20,25 Z" />
        </svg>
        ''',
        'style1': '''
        <svg viewBox="0 0 50 50">
            <path fill="black" d="
                M25,5 L40,25 L25,45 L10,25 Z
                M25,15 L35,25 L25,35 L15,25 Z" />
        </svg>
        ''',
        'style2': '''
        <svg viewBox="0 0 50 50">
            <path fill="black" d="
                M25,10 L45,25 L25,40 L5,25 Z" />
        </svg>
        '''
    }
    return icons.get(style, icons['default'])

# Plot and save top 10 species
plt.figure(figsize=(12, 6))
bars = plt.barh(species_data['Taxa'], species_data['Larvae_Detections'], color='skyblue')
plt.xlabel('Number of Larvae')
plt.ylabel('Species')
plt.title('Top 10 Most Detected Species in Guma\n(excluding nests S & T)')
plt.gca().invert_yaxis()

# Add value labels on bars
for bar in bars:
    width = bar.get_width()
    plt.text(width, bar.get_y() + bar.get_height()/2, 
             f'{int(width)}', 
             ha='left', va='center', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{save_dir}/top_10_species_larvae_counts.png', 
            dpi=300, bbox_inches='tight')
plt.show()

# Plot and save top 10 genera with icons
plt.figure(figsize=(12, 6))
bars = plt.barh(genera_data.index, genera_data.values, color='salmon')
plt.xlabel('Number of Larvae')
plt.ylabel('Genus')
plt.title('Top 10 Most Detected Genera in Guma\n(excluding nests S & T)')
plt.gca().invert_yaxis()

# Add icons and numbers
for i, bar in enumerate(bars):
    width = bar.get_width()
    y_pos = bar.get_y() + bar.get_height()/2
    
    # Create icon with different style based on index
    style = f'style{(i % 2) + 1}' if i < 2 else 'default'
    icon = create_moth_icon(style)
    
    # Add icon annotation
    ax_icon = plt.gca().add_artist(
        plt.Circle((width + 2, y_pos), 0.3, 
                  facecolor='black', 
                  edgecolor='none'))
    
    # Add number below icon
    plt.text(width + 2, y_pos - 0.4, 
             f'{int(width)}', 
             ha='center', va='top', 
             fontsize=8)

plt.tight_layout()
plt.savefig(f'{save_dir}/top_10_genera_larvae_counts_with_icons.png', 
            dpi=300, bbox_inches='tight')
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print("-" * 50)
print(f"Total number of larvae analyzed: {len(guma_columns)}")
print(f"Total number of unique species: {sum(guma_data['Is_Species'])}")
print(f"Total number of unique genera: {len(genus_detections)}")

print("\nTop 10 Species Detections:")
print("-" * 50)
for _, row in species_data.iterrows():
    print(f"{row['Taxa']}: detected in {int(row['Larvae_Detections'])} larvae")

print("\nTop 10 Genera Detections:")
print("-" * 50)
for genus, count in genera_data.items():
    print(f"{genus}: detected in {int(count)} larvae")

# Save detection data to Excel
# Create DataFrames for export
species_export = pd.DataFrame({
    'Species': species_data['Taxa'],
    'Number_of_Larvae': species_data['Larvae_Detections']
})

genera_export = pd.DataFrame({
    'Genus': genera_data.index,
    'Number_of_Larvae': genera_data.values
})

# Save to Excel
with pd.ExcelWriter(f'{save_dir}/detection_counts.xlsx') as writer:
    species_export.to_excel(writer, sheet_name='Top_10_Species', index=False)
    genera_export.to_excel(writer, sheet_name='Top_10_Genera', index=False)