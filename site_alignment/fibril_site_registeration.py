# This script is to run the fibril site registeration 

## libraries
import os, glob, argparse, datetime
import numpy as np
import pandas as pd
from fibrilsite import execute_global_registration, refine_registration, ply_parser_hull, best_registration, global_reg_pipeline, registrate_all_pockets

## Functions
def create_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parse = argparse.ArgumentParser()
    parse.add_argument("--info_file",     type=str, nargs=1, help="CSV containing the sites input feats.")
    parse.add_argument("--sites_folder",  type=str, nargs=1, help="Folder contaning site files.")
    parse.add_argument("--output_folder", type=str, nargs=1, help="Path to results export folder.")
    return parse

def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args

## Execution
def main():
    
    """ Main execution point """
    
    # Parse arguments
    args = parse_args(create_parser())
    df_pockets  = pd.read_csv(os.path.abspath(args.info_file[0]), index_col=0).rename(columns={'pocket_id':'pocket'})
    files       = os.path.abspath(args.sites_folder[0])
    output_root = os.path.abspath(args.output_folder[0])

    ## make output dir
    output = os.path.join(os.path.abspath(output_root), str(datetime.date.today())+"_registration_outputs")
    os.makedirs(output, exist_ok=0)
    
    # check point
    print(f"loaded df_pockets {df_pockets.shape}")
    
    ### get the path to ply files of each pocket
    path_dic = {}
    
    for path in glob.iglob(os.path.join(files, "*", "*_hull.ply")):
        path_dic[os.path.basename(path).split('convex')[0].strip('_')] = path
    
    ### creating a list of the pocket names based on the directory dictionnary
    names = np.array(list(path_dic.keys()))
    #print(names)
        
    #### perform alignments
    """ 
    The pipeline has the following characteristics:
        - Uses input features (including shape index) as features
        - Matches on the correspondence set
    """
    df_input = registrate_all_pockets(n_regs=5, path_dic=path_dic, df_pockets=df_pockets, output=output)
    df_input.to_csv(f'{output}/{datetime.date.today()}_all_pockets_input_features_registered.csv')
    
    print('Alignments done !')
    
    return None

if __name__ == "__main__":
    main()