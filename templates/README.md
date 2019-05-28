# Template Extractor
--------------------

Code to clean and extract templates from USPTO reaction database.

### Getting Started

##### Pre-requisites

* 7z archive extraction tool
* rdkit
* numpy
* pandas
* joblib


##### Step 1

Download USPTO reaction database from https://figshare.com/articles/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873, extract 7z archive, and place `1976_Sep2016_USPTOgrants_smiles.rsmi` into `data/` folder

```bash
$ mkdir data/ && cd data/
$ wget -O 1976_Sep2016_USPTOgrants_smiles.7z https://ndownloader.figshare.com/files/8664379
$ 7z e 1976_Sep2016_USPTOgrants_smiles.7z
$ cd ../
```

##### Step 2

Run `clean_and_extract_uspto.py` script. This will try to use all the CPU cores on your machine. On 32 cores it takes roughly 1 hour to run.

```bash
$ python clean_and_extract_uspto.py
```

This will generate `data/uspto.reactions.json.gz` and `data/uspto.templates.json.gz`. These two files can also be downloaded directly from [here](https://chemrxiv.org/articles/RDChiral_An_RDKit_Wrapper_for_Handling_Stereochemistry_in_Retrosynthetic_Template_Extraction_and_Application/7949024) if you do not wish to re-run the extraction code.

TODO:
* better documentation
* argparse
* template grouping
* merge templates with reactions to generate training data
