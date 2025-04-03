# Libraries
import os, glob, math, statistics, shutil, copy, prody, datetime

import pymesh
import open3d as o3d

import numpy as np
import pandas as pd
from shutil import copyfile

import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

from SBI.structure import PDB, Frame3D

from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

# Functions
## File parsers
def pdb_parser_and_per_atom_sasa_calculator(out_path, pdb_file, probe_r=1.4,):
    """
    Function will parse the pdb file, calculate the sasa and return a dataframe with atomid, coords and sasa
    
    out_path : File exporting desitination
    pdb_file : pdb file to be parsed
    probe_r  : The radius of the probe used for sasa calculation
    """
     
    # info dict
    coords_sasa = {} 
    coords_sasa['chain']     = []
    coords_sasa['resid']     = []
    coords_sasa['resname']   = []
    coords_sasa['atom_type'] = []
    coords_sasa['coords']    = []
    coords_sasa['sasa']      = []
    
    # parse pdb file
    parser = PDBParser(QUIET=1)
    struct = parser.get_structure(os.path.basename(pdb_file).split('_')[0], pdb_file)[0]
    chains = [chain.id for chain in struct]
    
    # calculate per atom sasa
    r = probe_r
    sr = ShrakeRupley(probe_radius=r, n_points=100,)
    sr.compute(struct, level='A')
    
    # Parse info
    for chain in chains:
        residues = [res.id[1] for res in struct[chain].get_residues()]
        for resid in residues:
            atoms = [atom.get_full_id() for atom in struct[chain][resid].get_atoms()]
            for atm in atoms:
                coords_sasa['chain'].append(atm[2])
                coords_sasa['resid'].append(atm[3][1])
                coords_sasa['resname'].append(struct[chain][resid].get_resname())
                coords_sasa['atom_type'].append(atm[4][0])
                coords_sasa['coords'].append((struct[chain][resid][atm[4][0]].get_coord()))
                coords_sasa['sasa'].append(round(struct[chain][resid][atm[4][0]].sasa,2))
    
    # create the coords and sasa data frame (for application)
    df_coords_sasa = pd.DataFrame(coords_sasa)
    df_coords_sasa.to_csv(f"{out_path}/{datetime.date.today()}_{os.path.basename(pdb_file).split('_')[0]}_per_atom_sasa_coords.csv")
    
    return df_coords_sasa

def ply_parser(out_path, ply_file, ):
    """
    This function is to parse MaSIF ply files and return a dataframe with the surface point coords, surface points normals and averaged shape index
    
    out_path : File exporting desitination
    ply_file : ply file to be parsed
    """
    # info container
    dic = {'surf_coords' :[],
           'surf_charge' :[],
           'surf_hbond'  :[],
           'surf_hphob'  :[],
           'surf_normals':[],
          }
    # Parse the ply file
    with open(ply_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if not line[0].isalpha() and len(line.split()) == 10]
    
    for line in lines:
        # get surface point coords
        surf_x, surf_y, surf_z = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
        surf_coords = np.array([surf_x, surf_y, surf_z])
        
        # get coord charge
        surf_charge = float(line.split()[4])
        
        # get coord hbond
        surf_hbond = float(line.split()[5])
        
        # get coord hphob
        surf_hphob = float(line.split()[6])
        
        # get surface normals
        surf_nx, surf_ny, surf_nz = float(line.split()[7]), float(line.split()[8]), float(line.split()[9])
        surf_norm_coords = np.array([surf_nx, surf_ny, surf_nz])
                
        dic['surf_coords'].append(surf_coords)
        dic['surf_charge'].append(surf_charge)
        dic['surf_hbond'].append(surf_hbond)
        dic['surf_hphob'].append(surf_hphob)
        dic['surf_normals'].append(surf_norm_coords)
    
    df_surf = pd.DataFrame.from_dict(dic)
    df_surf.reset_index(drop=False, inplace=True)
    df_surf.rename(columns={'index':'MaSIF_index'}, inplace=True)
    df_surf.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(ply_file).split("_")[0]}_masif_params_from_ply.csv')
    
    return df_surf

## Vector calculators
def calculate_eigenvectors(pdb_file, anchor_residues, anchor_chains):
    """
    This function is to calculate fibril eigenvectors
    
    pdb_file        : file to be used
    anchor_residues : the residue(s) index that will be used to calculate the fibril eigenvectors, MUST BE A LIST
    anchor_chains   : the chains to be used for calculating the fibril eigenvectors, MUST BE A LIST
    """
    # Load and clean the pdb file
    pdb3d  = PDB(pdb_file, format='pdb', clean=True,
                     dehydrate=True, hetatms=False)#['AtomTask:PROTEINBACKBONE']
    
    # Get the fibril part
    chunk = pdb3d[(pdb3d.auth_seq_id.isin(anchor_residues)) & (pdb3d.auth_asym_id.isin(anchor_chains))]
    
    # Calculate eigenvectors
    evecs = chunk.eigenvectors(10)
    side, perpend, major = evecs[0][-1] - evecs[0][0], evecs[1][-1] - evecs[1][0], evecs[2][-1] - evecs[2][0]
    side    = side    / np.linalg.norm(side)
    perpend = perpend / np.linalg.norm(perpend)
    major   = major   / np.linalg.norm(major)
    
    return side, perpend, major

## Isolate the fibril side surface points
def isolate_fibril_side(out_path, pdb_file, df, fibril_axis):
    """
    This function is to get the surface points on the sides of the fibrils only using normals calculations
    
    out_path    : File exporting desitination
    pdb_file    : pdb file for the fibril
    df          : Dataframe containing the information
    fibril_axis : Fibril elongation axis (major)
    """
    # get dot product of the fibril axis with the surface nomals
    df['surf_norm_fibril_dot'] = df['surf_normals'].apply(lambda x: np.dot(fibril_axis,x))
    
    # exclude top and bottom surface points and get fibril side surface only, these values were set based on experimentation
    df_side_surf = df[(df['surf_norm_fibril_dot']  <= 0.3) & (df['surf_norm_fibril_dot']  >= -0.3)]
    df_side_surf.reset_index(drop=True, inplace=True)
    
    # export
    df_side_surf.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(pdb_file).split("_")[0]}_fibril_side_surface_points.csv')
    
    return df_side_surf

## Isolate groove surface points from atomic coordinates
def pocket_getter(df, pocket_name, pdb_file, out_path, sym=True, chains_sym=None, chains_sym_resids=None,
                  chains_asym_1=None, chains_asym_1_resids=None,
                  chains_asym_2=None, chains_asym_2_resids=None, ):
    """
    This function is to get the properties for a pocket given its forming residues
    
    df                   : Dataframe containing information
    pocket_name          : Pocket identifier
    pdb_file             : pdb file for the fibril
    out_path             : File exporting desitination
    sym                  : True --> if the groove is in one protofilament and False --> if the spanned within different protofilament
    chains_sym           : Chains composing the groove
    chains_sym_resids    : Resids in chains_sym
    chains_asym_1        : Chains comoposing a part of the groove
    chains_asym_1_resids : Resids in chains_asym_1 
    chains_asym_2        : Chains comoposing a part of the groove
    chains_asym_2_resids : Resids in chains_asym_2    
    """
    
    dic = {}                  # dic for composing the chain-resid pair
    res_chain_pairs = []      # pair res to chains
    vessel = []               # container for filtered df
    
    if sym == True and (chains_sym != None and chains_sym_resids != None):
        df_sel_res = df[(df.resid.isin(chains_sym_resids)) & (df.chain.isin(chains_sym))]
        df_sel_res.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(pdb_file).split("_")[0]}_{pocket_name}_coords.csv')
    
        return df_sel_res
    
    elif sym == False and (chains_asym_1 != None and chains_asym_1_resids != None) and (chains_asym_2 != None and chains_asym_2_resids != None):
        for chain in chains_asym_1:
            dic[chain] = chains_asym_1_resids
        for chain in chains_asym_2:
            dic[chain] = chains_asym_2_resids
        for c in dic.keys():
            for r in dic[c]:
                res_chain_pairs.append((c,r))
        for i, chain, resid in zip(df.index, df.chain, df.resid):
            if (chain, resid) in res_chain_pairs:
                vessel.append(df.iloc[i])
        df_pocket_coords = pd.DataFrame(vessel)
        df_pocket_coords.reset_index(drop=True, inplace=True)
        # export
        df_pocket_coords.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(pdb_file).split("_")[0]}_{pocket_name}_coords.csv')
        
        return df_pocket_coords

## Surface to Coords mappers and etc.
def surface_atomic_coords_mapper(out_path, df_surf, df_atomic, pdb_file, threshold=5.0,):
    """
    This function is to map the surface points to the closest atomic coordinates
    Surface points coming from the MaSIF ply file --> Shape index, patch retrival
    Atmoic coordinates coming from the pdb file --> per atom sasa
    
    out_path    : File exporting desitination
    df_surf     : Dataframe containing information of surface points
    df_atomic   : Dataframe containing information of atomic coordinates
    threshold   : Threshold for connecting the surface and atomic coords
    pdb_file    : pdb file for the fibril
    """
    
    # collection container
    vessel = []
    
    #getting the coords
    surf_points = np.array([item for item in df_surf['surf_coords']])
    atm_coords = np.array([item for item in df_atomic['coords']])
    
    #generate the distance matrix
    dists = cdist(atm_coords,surf_points)
    
    #get the points that are closest to each other 
    closest_points = np.argmin(dists, axis=1)
    
    # filter by threshold
    closest_coords = closest_points[dists[np.arange(len(dists)), closest_points] < threshold]
    
    # get the atom indicies closest to the surface points as they are arranged according to the surface points indicies
    atom_idx = closest_coords.tolist()
    
    for points, close_idx in zip(range(len(df_atomic.index)), atom_idx):
        vessel.append(df_surf.iloc[close_idx])

    df_closest_pc = pd.DataFrame(vessel)
    df_closest_pc.reset_index(drop=True, inplace=True)

    # The merge
    df_merged = pd.merge(df_atomic, df_closest_pc, how='inner', left_index=True, right_index=True)
    df_merged.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(pdb_file).split("_")[0]}_mapped_coords_surf_and_atomic.csv')
    
    return df_merged

def surface_atomic_coord_finder(df_surf, atomic_coord, pdb_file,):
    """
    This function is to map the surface points to the closest atomic coordinates
    Surface points coming from the MaSIF ply file --> Shape index, patch retrival
    Atmoic coordinates for certain atom
    
    out_path    : File exporting desitination
    df_surf     : Dataframe containing information of surface points
    atomic      : Atomic coord to find the closest surface point for
    pdb_file    : pdb file for the fibril
    """
    
    # collection container
    vessel = []
    
    #getting the coords
    surf_points = np.array([item for item in df_surf['surf_coords']])
    atm_coords = atomic_coord[None,:]
    
    #generate the distance matrix
    dists = cdist(atm_coords,surf_points)
    
    #get the points that are closest to each other 
    closest_points = np.argmin(dists, axis=1)
    
    # get the atom indicies closest to the surface points as they are arranged according to the surface points indicies
    atom_idx = closest_points.tolist()
    
    return df_surf.iloc[atom_idx]

def surface_coord_close_neighbours(df_surf, surf_coord, pdb_file, current_points, threshold=5.0,):
    """
    This function is to map the surface points to the closest atomic coordinates
    Surface points coming from the MaSIF ply file --> Shape index, patch retrival
    Atmoic coordinates for certain atom
    
    out_path       : File exporting desitination
    df_surf        : Dataframe containing information of surface points
    surf_coord     : Surface coord to find the closest surface point for
    pdb_file       : pdb file for the fibril
    current_points : list of current surf points for a check to avoid redunduncy
    threshold      : Threshold for connecting the surface and atomic coords in Å
    """
    
    # collection container
    vessel = []
    
    #getting the coords
    surf_points = np.array([item for item in df_surf['surf_coords']])
    surf_coords = surf_coord[None, :]
    
    #generate the distance matrix
    dists = cdist(surf_coords, surf_points)
    
    closest_points = np.argwhere(dists < threshold) 
    
    ## get the atom indicies closest to the surface points as they are arranged according to the surface points indicies
    surf_idx = [idx_df for idx_coord,idx_df in closest_points]
    
    for idx in surf_idx:
        if list(df_surf.iloc[idx].surf_coords) not in current_points:
            vessel.append(df_surf.iloc[idx])
    
    df = pd.DataFrame(vessel)
        
    return df

## The Groove Fxs

def the_groove_expander(df_pocket_surf, df_atom_sasa, df_side_surf, pdb_file, out_path, pocket_name, threshold=3.0,):
    """
    This function is add all the point forming the pocket that were not captured during the groove walking 

    df_pocket_surf      : Dataframe containing the mapped coords to surface points
    df_atom_sasa        : Dataframe contianing the atoms' coords
    df_side_surf        : Dataframe containing the isolated fibril side surface point
    pdb_file            : pdb file path
    out_path            : File exporting desitination
    pocket_name         : Pocket identifier
    threshold           : Threshold to be used for calculating the closest neighbours to a point, default 3Å 
    """ 
    # Get a set of surf_coords
    isolated_surf_coords = [list(coord) for coord in df_pocket_surf.surf_coords]
    
    # Containers
    to_keep = []
    vessel = []
    
    # Expanding coords isolation
    for surf_coord, surf_norm, current_chain in zip (df_pocket_surf.surf_coords, df_pocket_surf.surf_normals, df_pocket_surf.chain):
        
        df_temp = surface_coord_close_neighbours(
            df_surf=df_side_surf,
            surf_coord=surf_coord,
            pdb_file=pdb_file,
            current_points=isolated_surf_coords,
            threshold=threshold,)

        vessel.append(df_temp)
        df_to_add = pd.concat(vessel)
        df_to_add['point_direction'] = df_to_add['surf_normals'].apply(lambda x: np.dot(x, surf_norm))
        df_to_keep = df_to_add[df_to_add.point_direction > 0.000]
        to_keep.append(df_to_keep)
    
    df_pocket_extension = pd.concat(to_keep)
    df_pocket_extension.reset_index(drop=False, inplace=True)
    df_pocket_extension.drop_duplicates(subset='index', inplace=True)
    df_pocket_extension.drop(columns='index',  inplace=True)
    df_pocket_extension.reset_index(drop=True, inplace=True)
    
    # Add to points to the pocket
    df_full_pocket = df_pocket_surf.append(df_pocket_extension)
    df_full_pocket.reset_index(drop=True, inplace=True)
    
    # export
    df_full_pocket.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(pdb_file).split("_")[0]}_{pocket_name}_pocket_after_expansion.csv')
    
    print(df_full_pocket.shape)
    
    return df_full_pocket

def pocket_rim_cleaner(df_atom_sasa, df_full_pocket, pocket_chains, pdb_file,
                             res_resid, res_atom, guard_resid, guard_atom, outpath, threshold =10.0):
    """
    This function is to clean up the pocket rims
    
    df_atom_sasa         : Dataframe contianing the atoms' coords
    df_full_pocket       : Dataframe containing the expanded groove coords
    pocket_chains        : Chain ids for the end of the pocket
    pdb_file             : pdb file path
    res_resid            : The resid for the first groove res
    res_atom             : The selected atom from the first res 
    guard_resid          : The residue outside the pocket to use for removing outside surface points
    guard_atom           : The defined atom for the guard_start_resid
    out_path             : File exporting desitination
    threshold            : Distance threshold from the guard residue to remove the outside points
    """
    
    # containers
    vessel = []

    # isolate the start rim surf points
    for chain in pocket_chains:
        # get the start residue surf coord
        atom_coord  = np.array(df_atom_sasa[(df_atom_sasa.resid == res_resid)& (df_atom_sasa.chain == chain) & (df_atom_sasa.atom_type == res_atom)]['coords'].tolist())[0]
        surf_coord  = np.array(surface_atomic_coord_finder(df_surf=df_full_pocket, atomic_coord=atom_coord, pdb_file=pdb_file)['surf_coords'].tolist())[0]
        
        # get surf coords for the rim guard 
        rim_guard_coord = np.array(df_atom_sasa[(df_atom_sasa.resid == guard_resid) & (df_atom_sasa.chain == chain) & (df_atom_sasa.atom_type == guard_atom)]['coords'].tolist())[0]
        
        df_test = df_full_pocket.copy()
        
        df_test['dist_from_guard'] = df_test['surf_coords'].apply(lambda x: np.linalg.norm(x-rim_guard_coord))
        df_test['dist_from_start'] = df_test['surf_coords'].apply(lambda x: np.linalg.norm(x-surf_coord))
        
        cleaned = df_test[df_test.dist_from_guard < threshold]
        vessel.append(cleaned)
        
    df = pd.concat(vessel)
    print(df.shape)
    
    return df

def get_areas(df):
    """
    This function is to calculate the protein surface area based on the per atom SASA. Also, this will calculate the hydrophobic portion of that surface
    df : Data frame containing SASA information
    """
    
    #defining the hydrophobic residues
    hphob_res  = ['ALA', 'ILE', 'LEU', 'MET', 'VAL', 'PHE', 'TRP', 'TYR', 'GLY',  'PRO', 'CYS']
    
    # Get the surface area estimate from calculated SASA
    total_sasa = round(df.sasa.sum(), 2)
    
    # Get the hydrophobic fraction SASA
    pocket_hphob_sasa = 0

    pocket_resnames = list(set(df.resname))     # get the pocket residue set
    grouped = df.groupby(by='resname')          # group by the resname
    
    for item in pocket_resnames:
        if item in hphob_res:
            df_temp = grouped.get_group(item)
            pocket_hphob_sasa += df_temp.sasa.sum()
            
    pocket_hphob_sasa = round(pocket_hphob_sasa, 2)
    percent_hphob = round((pocket_hphob_sasa/total_sasa)*100, 2)
    
    return total_sasa, pocket_hphob_sasa, percent_hphob

## PyMOL visualization
def pymol_xyz_writer(out_path, file_name, df, coords_column, element='C'):
    """
    This function is to write a pymol xyz file for points visualization
    
    out_path      : File exporting desitination
    file_name     : File name
    df            : Dataframe containing information
    coords_column : Name of column containing x,y,z coordinates
    element       : The element representing the surface points in question
    """
    with open(f"{out_path}/{datetime.date.today()}_{file_name}.xyz", "w") as w:
        for coords in df[f'{coords_column}']:
            w.write(f' {element} {coords[0]} {coords[1]} {coords[2]} \n')
    
    return 'GO SEE!'

def vector_xyz_writer(out_path, file_name, vector, element='C'):
    """
    This function is to write a pymol xyz file for points visualization
    
    out_path      : File exporting desitination
    file_name     : File name
    vector        : Axis x,y,z coordinates
    element       : The element representing the surface points in question
    """
    
    eigenvec = vector[None,:] * np.linspace(0,1,50)[:,None] * 35
    
    with open(f"{out_path}/{datetime.date.today()}_{file_name}.xyz", "w") as w:
        for coords in eigenvec:
            w.write(f' {element} {coords[0]} {coords[1]} {coords[2]} \n')
    
    return 'GO SEE!'

## Additional manual filtering
def filter_by_range(df_atom_sasa, df_pocket_isolate, resid, chain, atom, output, pocket_id, pdb_file, dist_thresh=5.0):
    """
    This function is to remove points in the proximity of a certain residue
    
    df_atom_sasa         : Dataframe contianing the atoms' coords
    df_pocket_isolate    : Dataframe containing the isolated surface point of the pocket
    resid                : Desired res as anchor
    chain                : Desired res's chain
    atom                 : Desired res's atom
    output               : File exporting desitination
    pocket_id            : Pocket name
    pdb_file             : pdb file path
    dist_thresh          : Distance from the res to remove the surface point
    """
    
    atom_coord  = np.array(df_atom_sasa[(df_atom_sasa.resid == resid) & (df_atom_sasa.chain == chain) & (df_atom_sasa.atom_type == atom)]['coords'].tolist())[0]
    surf_coord  = np.array(surface_atomic_coord_finder(df_surf=df_pocket_isolate, atomic_coord=atom_coord, pdb_file=pdb_file)['surf_coords'].tolist())[0]
    surf_norm   = np.array(surface_atomic_coord_finder(df_surf=df_pocket_isolate, atomic_coord=atom_coord, pdb_file=pdb_file)['surf_normals'].tolist())[0]
        
    df_to_remove = df_pocket_isolate.copy()
    
    df_to_remove['dist_from_mid']            = df_to_remove['surf_coords'].apply(lambda x: np.linalg.norm(x-surf_coord))
    df_to_remove['point_norm_dot_mid_norm']  = df_to_remove['surf_coords'].apply(lambda x: np.dot(x,surf_norm))
    
    df_to_remove = df_to_remove[(df_to_remove.dist_from_mid < dist_thresh)]
    
    print(df_to_remove.shape)
    idx = df_to_remove.index.tolist()
    
    df_pocket_refined = df_pocket_isolate.copy()
    df_pocket_refined = df_pocket_refined[~df_pocket_refined.index.isin(idx)]
    df_pocket_refined.reset_index(drop=True, inplace=True)
    
    print(df_pocket_refined.shape)
    
    # export
    df_pocket_refined.to_csv(f'{output}/{datetime.date.today()}_{os.path.basename(pdb_file).split("_")[0]}_{pocket_id}_refined.csv')
    
    return df_pocket_refined

# Alignmnet functions
def execute_global_registration(source_pcd, target_pcd, source_feats,
                                target_feats, distance_threshold = 1.5, verbose:bool = False):
    """
    This function is to use Ransac to align the pockets' point clouds
    
    source_pcd           : source pcd that will be the target to align on
    target_pcd           : target pcd that will be aligned to source pcd
    source_feats         : source pcd features that will be used in the alignment
    target_feats         : target pcd features that will be used in the alignment
    distance threshold   : Ransac distance threshold to align two pcd, 1.5 was adopted from MaSIF seed search  
    verbose                    : Parameter controlling verbosity
    """
    
    if verbose:
        print(":: RANSAC registration on downsampled point clouds.")
        print("   we use a liberal distance threshold %.3f." % distance_threshold)
    
    result = o3d.pipelines.registration.registration_ransac_based_on_feature_matching(
        source_pcd, target_pcd, source_feats, target_feats, True,
        distance_threshold,
        o3d.pipelines.registration.TransformationEstimationPointToPoint(False),
        3, 
        [
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnEdgeLength(0.9),
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnDistance(distance_threshold),
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnNormal(np.pi/2),
        ], 
        o3d.pipelines.registration.RANSACConvergenceCriteria(20000, 0.999))
    return result

def refine_registration(source_pcd, target_pcd, result_ransac, distance_threshold = 1.0, verbose:bool = False, generalized: bool = False):
    
    """
    This function is to use ICP to refine the aligned the pockets' point clouds
    
    source_pcd           : source pcd that will be aligned to the target pcd (apply the transformation matrix to this pcd)
    target_pcd           : target pcd that the source pcd will aligned to
    source_feats         : source pcd features that will be used in the alignment
    target_feats         : target pcd features that will be used in the alignment
    distance threshold   : Ransac distance threshold to align two pcd, 1.0 was adopted from MaSIF seed search  
    """
    
    if verbose:
        print(":: Point-to-plane ICP registration is applied on original point")
        print("   clouds to refine the alignment. This time we use a strict")
        print("   distance threshold %.3f." % distance_threshold)
    
    if generalized:
        result = o3d.pipelines.registration.registration_generalized_icp(
            source_pcd, target_pcd, distance_threshold, result_ransac.transformation,
        )
    
    else: 
        result = o3d.pipelines.registration.registration_icp(
            source_pcd, target_pcd, distance_threshold, result_ransac.transformation,
            o3d.pipelines.registration.TransformationEstimationPointToPlane(),
            o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=100)
        )
    return result

def ply_parser_hull(ply_file:str, out_path=None ) -> pd.DataFrame:
    """
    This function is to parse convex hull MaSIF ply files and return a dataframe with the surface point coords and surface points normals.
    
    out_path : File exporting desitination
    ply_file : ply file to be parsed
    """
    # info container
    dic = {'surf_coords' :[],
           'surf_normals':[],
          }
    # Parse the ply file
    with open(ply_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if not line[0].isalpha() and len(line.split()) == 10]
    
    for line in lines:
        # get surface point coords
        surf_x, surf_y, surf_z = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
        surf_coords = np.array([surf_x, surf_y, surf_z])
                
        # get surface normals
        surf_nx, surf_ny, surf_nz = float(line.split()[7]), float(line.split()[8]), float(line.split()[9])
        surf_norm_coords = np.array([surf_nx, surf_ny, surf_nz])

        # add to dict       
        dic['surf_coords'].append(surf_coords)
        dic['surf_normals'].append(surf_norm_coords)
    
    df_surf = pd.DataFrame.from_dict(dic)
    df_surf.reset_index(drop=False, inplace=True)
    df_surf.rename(columns={'index':'MaSIF_index'}, inplace=True)
    if out_path != None:
        df_surf.to_csv(f'{out_path}/{datetime.date.today()}_{os.path.basename(ply_file).split("_")[0]}_masif_params_from_ply.csv')
    
    return df_surf

def best_registration(source_pcd:o3d.geometry.PointCloud, target_pcd:o3d.geometry.PointCloud, source_feats:o3d.pipelines.registration.Feature, target_feats:o3d.pipelines.registration.Feature, n:int = 3, verbose:bool = False):
    """Takes the source and target point cloud, and performs n instances of global registration. In cases where the n 
    attempts fail to find any correspondence, the process is repeated. Outputs the registration result that yields the biggest
    number of correspondences
    
        Args:
        source_pcd (o3d.geometry.PointCloud): Source point cloud
        target_pcd (o3d.geometry.PointCloud): Target point cloud
        source_feats (o3d.pipelines.registration.Feature): Source features
        target_feats (o3d.pipelines.registration.Feature): Target Features
        n (int, optional): Number of global registrations to perform. Defaults to 3.
        verbose (bool, optional): Whether or not to print the performance of the registration. Defaults to False

    Returns:
        (o3d.pypeline.registration.registration_results): _description_ [0] for ransac results, [1] for icp results
    """
    # containers
    results_ransac  = []
    results_icp     = []
    
    n_corres_ransac = []
    n_corres_icp    = []
    
    measures_icp    = []

    n_max = 0

    while n_max == 0:
        for i in range(n):
            result_ransac = execute_global_registration(
                source_pcd=source_pcd,
                target_pcd=target_pcd,
                source_feats=source_feats,
                target_feats=target_feats,
                distance_threshold= 1.5,
            )
            result_icp = refine_registration(
                source_pcd=source_pcd,
                target_pcd=target_pcd,
                result_ransac=result_ransac,
                distance_threshold=1
            )
            
            results_ransac.append(result_ransac)
            results_icp.append(result_icp)
            
            n_corres_ransac.append(len(np.asarray(result_ransac.correspondence_set)))
            n_corres_icp.append(len(np.asarray(result_icp.correspondence_set)))
            
            measures_icp.append([result_icp.inlier_rmse, result_icp.fitness])   # for the verbose function
            
        n_max = max(n_corres_icp)
    
    if verbose:
        print('Chosen alignment measures', measures_icp[np.argmax(n_corres_icp)])
        
    return results_ransac[np.argmax(n_corres_icp)], results_icp[np.argmax(n_corres_icp)]

def global_reg_pipeline(source_pocket:str, target_pocket:str, path_dic:dict, df_pockets:pd.DataFrame, output:str, n_regs:int = 5 ):
    """This function takes the name of two pockets (source and target) and performs global registration (RANSAC followed by ICP)
        It outputs different metrics on the registration

    Args:
        source_pocket (string): Name of the source pocket
        target_pocket (string): Name of the target pocket
        n_regs (int, optional): Number of time to perform the global registration for each pair of pockets. Defaults to 5
        matching_on_corres (bool, optional): What type of point matching to use. Defaults to True

    Returns:
        pd.DataFrame: Contains valuable metrics on the registration process

    Modifications:
        Removed "matching on all" doesn't make sense the way it was coded, you need a distance cutoff otherwise the matching can be very off, besides one target point can be matched to multiple source points.
    """

    # format paths
    source_path = f"{path_dic[source_pocket]}"
    target_path = f"{path_dic[target_pocket]}"

    # parse ply files
    df_source = ply_parser_hull(ply_file=source_path)
    df_target = ply_parser_hull(ply_file=target_path)

    target_surf_normals = df_target.surf_normals.to_numpy()
    source_surf_normals = df_source.surf_normals.to_numpy()

    # fetching the input feats from the pocket DataFrame
    source_input_si     = df_pockets[df_pockets['pocket'] == source_pocket]['input_si'].to_numpy()
    source_surf_charge  = df_pockets[df_pockets['pocket'] == source_pocket]['input_charge'].to_numpy()
    source_surf_hphob   = df_pockets[df_pockets['pocket'] == source_pocket]['input_hphob'].to_numpy()
    source_surf_hbonds  = df_pockets[df_pockets['pocket'] == source_pocket]['input_hbonds'].to_numpy()

    target_input_si     = df_pockets[df_pockets['pocket'] == target_pocket]['input_si'].to_numpy()       
    target_surf_charge  = df_pockets[df_pockets['pocket'] == target_pocket]['input_charge'].to_numpy()       
    target_surf_hphob   = df_pockets[df_pockets['pocket'] == target_pocket]['input_hphob'].to_numpy()       
    target_surf_hbonds  = df_pockets[df_pockets['pocket'] == target_pocket]['input_hbonds'].to_numpy()         

    feats_source = np.stack(( source_surf_charge, source_surf_hphob, source_surf_hbonds, source_input_si), axis=1)
    feats_target = np.stack(( target_surf_charge, target_surf_hphob, target_surf_hbonds, target_input_si), axis=1)

    # get the features for the alignment, not the calculated ones
    source_feats       = o3d.pipelines.registration.Feature()
    source_feats.data  = feats_source.T

    target_feats       = o3d.pipelines.registration.Feature()
    target_feats.data  = feats_target.T


    ### --- preprocessing --- 
    
    # substract the center of mass of each PC from itself
    source = o3d.io.read_point_cloud(source_path, format='ply')
    target = o3d.io.read_point_cloud(target_path, format='ply')

    # remove center of mass
    source_corr = o3d.geometry.PointCloud()
    source_corr.points  = o3d.utility.Vector3dVector(np.asarray(source.points) - np.mean(np.asarray(source.points), axis=0))
    source_corr.normals = o3d.utility.Vector3dVector(source_surf_normals)

    target_corr = o3d.geometry.PointCloud()
    target_corr.points  = o3d.utility.Vector3dVector(np.asarray(target.points) - np.mean(np.asarray(target.points), axis=0))
    target_corr.normals = o3d.utility.Vector3dVector(target_surf_normals)

    # storing the array form of the source and target
    source_points = np.asarray(source_corr.points)
    target_points = np.asarray(target_corr.points)
    size_source = len(source_points)
    size_target = len(target_points)

    ### --- running the registration ---

    result_ransac, result_icp = best_registration(source_pcd=source_corr,
            target_pcd=target_corr,
            source_feats=source_feats,
            target_feats=target_feats,
            n = n_regs
    )

    # storing the performances and attributes of the results
    ## from the ransac run
    ransac_rmse = result_ransac.inlier_rmse
    ransac_fitness = result_ransac.fitness

    ## from the icp run
    icp_rmse = result_icp.inlier_rmse
    icp_fitness = result_icp.fitness
    
    # export the o3d registration object
    dic = {}
    dic['ransac_test'] = {'corres_set':[np.asarray(result_ransac.correspondence_set)],
                          'transformation':[result_ransac.transformation],
                          'fitness': [result_ransac.fitness],
                          'rmse': [result_ransac.inlier_rmse]}

    dic['icp_test'] = {'corres_set':[np.asarray(result_icp.correspondence_set)],
                        'transformation':[result_icp.transformation],
                        'fitness': [result_icp.fitness],
                        'rmse': [result_icp.inlier_rmse]}
    
    ## a dir for the npy files
    npy_out = f'{output}/o3d_objects_npy/'
    if not os.path.exists(npy_out):
        os.makedirs(npy_out)
    
    npy_out_input_feats = f'{npy_out}/input_feats_reg/'
    if not os.path.exists(npy_out_input_feats):
        os.makedirs(npy_out_input_feats)
    
    ## save the objects
    np.save(f"{npy_out_input_feats}/{source_pocket}_{target_pocket}_input-feats.npy", dic, allow_pickle=True)


    # creating a dataframe that includes all informations we want to analyse
    output_df = pd.DataFrame({
        'source_pocket':[source_pocket], 
        'target_pocket':[target_pocket], 
        'ransac_rmse':[ransac_rmse], 
        'ransac_fitness':[ransac_fitness],
        'ransac_nb_corres':[len(np.asarray(result_ransac.correspondence_set))],
        'icp_rmse':[icp_rmse], 
        'icp_fitness':[icp_fitness], 
        'icp_nb_corres':[len(np.asarray(result_icp.correspondence_set))], 
        'size_source':[size_source], 
        'size_target':[size_target], 
    })

    return output_df

def registrate_all_pockets(n_regs:int, path_dic:dict, df_pockets:pd.DataFrame, output:str) -> pd.DataFrame:
    """Registrates all pockets

    Args:
        n_regs (int): See doc global_reg_pipeline()
    """
    # container for output df
    vessel = []

    # loop over every pocket in the training set (about 8min)
    for n in range(len(path_dic.keys())):
        for m in range(n+1, len(list(path_dic.keys()))):

            source = list(path_dic.keys())[n]
            target = list(path_dic.keys())[m]
            print(f"source:{source} , target:{target}")

            result = global_reg_pipeline(source_pocket=source, target_pocket=target, path_dic=path_dic, df_pockets=df_pockets, n_regs = n_regs, output=output)
            vessel.append(result)
    
    df = pd.concat(vessel).reset_index(drop=1)

    return df

# Alignmnet analysis functions

def get_input_desc_diff(df_all_info:pd.DataFrame, df:pd.DataFrame, idx:int, input_feats_npy:list) -> pd.DataFrame:
    """ This function calculates the Input features difference """ 
    
    # define the open3d results file path 
    path_ = [p for p in input_feats_npy if df.at[idx, 'source_pocket'] in p
            if df.at[idx, 'target_pocket'] in p]
    assert len(path_) == 1, 'Path finding ERROR'
    path_ = path_[0]
    
    # load the results file
    reg_results = np.load(path_, allow_pickle=True).item()['icp_test']
    
    # get the correspondences indicies
    source_corres_set = reg_results['corres_set'][0][:,0]
    target_corres_set = reg_results['corres_set'][0][:,1]
    assert len(source_corres_set) == len(target_corres_set)    
    
    # get the corres set properties
    ## for source pocket
    df_src = df_all_info[df_all_info.pocket_id == df.at[idx, 'source_pocket']].reset_index(drop=True)
    source_input_si    = df_src['input_si'].to_numpy()[source_corres_set]
    source_surf_charge = df_src['input_charge'].to_numpy()[source_corres_set]
    source_surf_hphob  = df_src['input_hphob'].to_numpy()[source_corres_set]
    source_surf_hbonds = df_src['input_hbonds'].to_numpy()[source_corres_set]
    
    feats_source = np.stack(( source_surf_charge, source_surf_hphob, source_surf_hbonds, source_input_si), axis=1)
    
    ## for target pocket
    df_target = df_all_info[df_all_info.pocket_id == df.at[idx, 'target_pocket']].reset_index(drop=True)
    target_input_si    = df_target['input_si'].to_numpy()[target_corres_set]
    target_surf_charge = df_target['input_charge'].to_numpy()[target_corres_set]
    target_surf_hphob  = df_target['input_hphob'].to_numpy()[target_corres_set]
    target_surf_hbonds = df_target['input_hbonds'].to_numpy()[target_corres_set]
    
    feats_target = np.stack(( target_surf_charge, target_surf_hphob, target_surf_hbonds, target_input_si), axis=1)
    
    assert feats_source.shape == feats_target.shape, 'Feats parsing failure'
    
    input_desc_diff = np.mean(np.linalg.norm(feats_target - feats_source, axis = 1))

    return input_desc_diff

def calc_input_feat_diff(df:pd.DataFrame, df_all_info:pd.DataFrame, input_feats_npy:list, output:str, export:bool=False, ident=None ) -> pd.DataFrame :
    """ This function is to calculate the input feat differece, raw and weighted """

    # calculate the input feature difference
    input_feat_diff_vessel = []

    for idx01 in df.index:
        input_feat_diff_v = get_input_desc_diff(df_all_info=df_all_info, df=df, idx=idx01, input_feats_npy=input_feats_npy)

        # append    
        input_feat_diff_vessel.append(input_feat_diff_v)

    # add to the dataframe
    df.insert(9, 'icp_mean_input_diff', input_feat_diff_vessel)

    assert len(input_feat_diff_vessel) == df.shape[0]

    if export and ident != None:
        df.to_csv(os.path.join(output, str(ident) + '_matches.csv'))
    
    return df