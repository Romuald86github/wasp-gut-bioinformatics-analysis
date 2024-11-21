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
   python tax_asv_merg.py


## 2. extracting the samples IDs from the cudadapt folder

9. run the script
    ```bash
    sample_IDs.py

## 3. obtain the taxonomy data mapped to their sample (gut) IDs 

10. run the script
    ```bash
    taxonomy_with_IDs.py



9. open the jupyter notebook interface by running ```jupyter notebook```
update 
copy the code in the script taxonomy_and_ASVs.py and change data paths to the relevant paths on your local machine
run the code in your jupyter notebook interface. 
