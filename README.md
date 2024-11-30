# wasp-gut-bioinformatics-analysis


# Bioinformatics-statistics-analysis

# I - data preparation

## 1 - merging taxonomic data to their original samples

after we screened out the wrongly detected species (not likely to be found in Sub-Saharan Africa), we then obtained the true taxonomic dataset, which was merged with ASV table containing ASV mapped with with the samples origins

### follow these steps to perform the merging:
open a terminal or command line interface

1. Clone the repository:
   ```bash
   git clone https://github.com/Romuald86github/Bioinformatics-statistics-analysis.git

conda create -n env

2. move into the repo

   ```bash
   cd wasp-gut-bioinformatics-analysis

   
2. create a conda virtual environment by running `(but first make sure you have anaconda installed on your computer)
   ```bash
   conda create -n env 


3. install pip with
   ```bash
   conda install pip
   
4. install jupyter notebook with
   ```bash
   pip install jupyter
   
6. install the python bioinformatics package
   ```bash
   pip install scikit-bio
   
8. move to data preparation
   ```bash
   cd 'data preparation'
9. run the script
   ```bash
   python tax_asv_merg.py \
   /path_to_your/table_ASVs.tsv \
   "path_to_your/Final taxonomy .csv" \
   /path_to_your/wasp_gut_results.tsv
  


## 2. extracting the samples IDs from the cudadapt folder

9. run the script
    ```bash
    python sample_IDs.py \
    /path/to/new/cutadapt/folder \
    /path/to/new/output/sample_origin_mapping_cleaned.tsv

## 3. obtain the taxonomy data mapped to their sample (gut) IDs 

10. run the script
    ```bash
    python taxonomy_with_IDs.py \
    /path/to/wasp_gut_results.tsv \
    /path/to/sample_origin_mapping_cleaned.tsv \
    /path/to/wasp_gut_results_with_origins.tsv


## 4. extracting metadata (larvae, nest, location)

11. run the script
    ```bash
    python metadata.py \
    /path/to/metadata.tsv \
    /path/to/wasp_gut_results_with_origins.tsv



## filtering out taxa detected with < 98% similraity

12. This was done using filter function in Excel (on wasp_gut_results_with_origins.tsv)


11. open the jupyter notebook interface by running () 
   ```bash
   jupyter notebook

 


# II - analysis (to be done within juoyter notebook)

1. plottinh accummulation curves


