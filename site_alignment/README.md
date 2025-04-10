# FibrilSite  – Site Alignment
To run the defined fibril sites alignment and alignment analysis:

## Provided scripts:
- 01_run_fibril_site_surface_feature_mapping.sh : \
    executable script to run [./scripts/fibril_site_input_feats_mapper.py] that would extract the surface features for the defined fibril site points following MaSIF processing of the fibril pdb files
    
    **params**: \
        - site_src_dir    : Path to folder containing the defined fibril sites. \
        - input_feats_src : Path to folder contaning MaSIF calculated surface features [/data_preparation/04b-precomputation_12A/precomputation] \
        - output_folder   : Path to root folder for exporting.

    **output**: \
        - {date}_all_sites_parsed.csv      : CSV containing information of all defined sites. \
        - {date}_all_sites_input_feats.csv : CSV containing information of all defined sites including the calculated surface features
    
    **exported information**: \
        - fibril               : fibril PDB id \
        - pocket_id            : pocket (site) name \
        - isolation            : level of pocket definition \
        - MaSIF_index          : pocket (site) surface points MaSIF index -> needed for mapping \
        - atom_type            : surface point atom type – if applicable \
        - chain                : fibril chain \
        - coords               : surface point atomic coordinates – if applicable \
        - point_direction      : dot product of surface normals of one anchor point and its neighbour -> needed for pocket expansion \
        - resid                : residue id \
        - resname              : residue name \
        - sasa                 : per residue solvent-accessible surface area \
        - surf_charge          : ply-file parsed surface electrostatics value \
        - surf_coords          : ply-file parsed surface coordinates \
        - surf_hbond           : ply-file parsed hydrogen-bond donor/acceptor potential \
        - surf_hphob           : ply-file parsed hydropathy value \
        - surf_norm_fibril_dot : dot product of point surface normal and fibril elongation axis \
        - surf_normals         : surface point normal \
        - input_si             : surface point calculated shape index \
        - input_charge         : surface point calculated electrostatics value \
        - input_hphob          : surface point calculated hydropathy value \
        - input_hbonds         : surface point calculated hydrogen-bond donor/acceptor potential 





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

