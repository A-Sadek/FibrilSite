# FibrilSite
**FibrilSite** is a computational pipeline to define fibril sites and identify shared features across a database of defined sites. 

Fibril surface sites are extracted from fibril surface as point clouds featurized with surface properties including Poisson–Boltzmann continuum electrostatics, hydrophobicity, shape index and hydrogen bond donors/acceptor patterns. 

### Fibril structures preparation
Fibril structures were pre-processed and surface features were computed using the **MaSIF** (Molecular Surface Interaction Fingerprinting) framework, as described in 
> Gainza, P. *et al.* *Nat Methods* **17**, 184–192 (2020). [https://doi.org/10.1038/s41592-019-0666-6]

### Fibril sites definition
Fibril site definition is peformed as described in the *fibril_grooves_extractor.ipynb* notebook in the example folder. 

To execute the notebook you need to install the code package containing all the required funtions as follows: \

1- create a conda environment with python version 3.6

    conda create -n fibrilsite python=3.6










### Fibril sites alignment 
