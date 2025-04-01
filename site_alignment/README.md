# FibrilSite  – Site Alignment

To run the defined sites alignment:

## Provided scripts:
 - fibril_site_registeration.py     : code for running the alignments
 - run_fibril_site_registeration.sh : script to run the alignments

## Usage

1- Create an inputs folder containing all the generated ply files for each site (e.g., Pol_1a_P1_convex_hull.ply)

2- Create a CSV file containing the information for the sites to be aligned using **fibril_site_input_feats_mapper.ipynb**

3- Adjust the parameters in the **run_fibril_site_registeration.sh** file

4- Activate the **fibrilsite** environment

5- Run the alignments as follows 
    
    source run_fibril_site_registeration.sh
