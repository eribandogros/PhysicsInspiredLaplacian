import numpy as np
import sys
import os
import argparse
import time
from biopandas.mol2 import PandasMol2
from joblib import Parallel, delayed

def get_xyz(mol2file):

    solvation = PandasMol2().read_mol2(mol2file).df

    solvation_e = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']

#    if elements == 'all':
#        solvation_e = ['H', 'C', 'N', 'O', 'F', 'S', 'Cl']
#    else:
#        solvation_e = re.findall('[A-Z][^A-Z]*',elements) #regex split by capital letter

    xyz_data = np.array([], dtype=float)
    for e in range(len(solvation_e)):
            ele = solvation_e[e]
            cloudpoint = \
                solvation[['x', 'y', 'z']][solvation['atom_name'].str.contains('^' + ele + '[0-9]*$')].values
            xyz_data = np.append(xyz_data, cloudpoint)

    return xyz_data

def convert_mol2(mol2_path,folder,id):
    mol2file=os.path.join(mol2_path, id + '.mol2')
    xyzfile=os.path.join(folder, id + '.npy')
    xyz_data=get_xyz(mol2file)
    np.save(xyzfile,xyz_data)
    
    
def main(args):
    dataset_name = args.dataset_name
    compound_list_file = args.compound_list_file
    prefix_name = args.prefix_name
    
    t1 = time.time()

    xyz_folder='../feature-gen/%s/xyz-folder-%s'%(dataset_name,prefix_name)
    if not os.path.exists(xyz_folder):
        os.mkdir(xyz_folder)

    dataset_location = '../'+dataset_name
    mol2_path = dataset_location+'/mol2-3D'
    list_compound_id=open(dataset_location+ '/'+ compound_list_file).read().splitlines()

    Parallel(n_jobs=30,backend="multiprocessing")(
        delayed(convert_mol2)(mol2_path,xyz_folder,id) for id in list_compound_id)

    print(f'Running time: {time.time() - t1}')


def parse_args(args):
    parser = argparse.ArgumentParser(description='Convert mol2 files in given dataset to xyz')

    parser.add_argument('--dataset_name', default='Dataset', type=str)
    parser.add_argument('--prefix_name', default='train', type=str)
    parser.add_argument('--compound_list_file', default='list_train.txt', type=str,
                        help='compound ids of a given dataset')
#    parser.add_argument('--elements', default='all', type=str,
#                        help='choose from H, C, N, O, F, S, Cl and concat or all')

    args = parser.parse_args()
    return args


def cli_main():
    args = parse_args(sys.argv[1:])
    print(args)
    main(args)


if __name__ == "__main__":
    cli_main()
    print('End!')
