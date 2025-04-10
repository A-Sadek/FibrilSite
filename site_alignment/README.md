# FibrilSite  – Site Alignment

To run the defined sites alignment:

## Provided scripts:
 - fibril_site_registeration.py     : code for running the alignments
 - run_fibril_site_registeration.sh : script to run the alignments
 - fibril_site_input_feats_mapper.ipynb : notebook to map the defined fibril site (i.e., pocket) points to their computed surface features

## Usage
1- Create a conda environment with python version 3.8
    
    conda create -n alignsites python=3.8 -y

2- Install the code package in the **fibrilsite** environment 

    pip install git+https://github.com/A-Sadek/FibrilSite.git

Or

    pip install fibrilsite

3- Install Openbabel

    conda install conda-forge::openbabel -y
    
4- Export the **alignsites** environment to Jupyter as follows:

    conda install -c anaconda ipykernel -y
    python -m ipykernel install --user --name=alignsites
 



1- Create an inputs folder containing all the generated ply files for each site (e.g., Pol_1a_P1_convex_hull.ply)

2- Create a CSV file containing the surface features information of the sites to be aligned using **fibril_site_input_feats_mapper.ipynb**

3- Adjust the parameters in the **run_fibril_site_registeration.sh** file

4- Activate the **alignsites** environment

5- Run the alignments as follows 
    
    source run_fibril_site_registeration.sh

