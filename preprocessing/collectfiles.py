import numpy as np
import sys
import os
import argparse

def check_for_missing_files(dataset_name,folder,dataset_id,prefix_name):
    missing_files = []
    for id in dataset_id:
         if not os.path.isfile(os.path.join(folder, id + '.txt')):
            missing_files.append(id.split('_').pop())
    if not missing_files == []:
        # save missing files to a .txt
        print(",".join(missing_files))
        with open("missing_files.txt", "w") as txt_file:
            txt_file.write(",".join(missing_files))
        print("Missing files, try again.")
        sys.exit()

def collect(dataset_name,folder,dataset_id,prefix_name,trim):
    features_collect = []
    n = 50
    for id in dataset_id:
        txt_feature_data = os.path.join(folder, id + '.txt')
        feature = np.loadtxt(txt_feature_data)
        if trim:
            # only want first n eigenvalues
            feature = feature[:n]
        features_collect.append(feature)
    features_collect = np.array(features_collect, dtype=object)
    
    if trim:
        np.save(f'{dataset_name}/{prefix_name}_{n}.npy',features_collect)
    else:
        np.save(f'{dataset_name}/{prefix_name}.npy',features_collect)

def main(args):
    dataset_name = args.dataset_name
    folder='../feature-gen/features/'
    trim = args.trim
    
    prefix_name = 'train'
    compound_list_file = 'list_'+prefix_name+'.txt'
    list_compound_id=open('../'+dataset_name+'/'+compound_list_file).read().splitlines()
    check_for_missing_files(dataset_name,folder,list_compound_id,prefix_name)
    collect(dataset_name,folder,list_compound_id,prefix_name,trim)
    
    prefix_name = 'test'
    compound_list_file = 'list_'+prefix_name+'.txt'
    list_compound_id=open('../'+dataset_name+'/'+compound_list_file).read().splitlines()
    check_for_missing_files(dataset_name,folder,list_compound_id,prefix_name)
    collect(dataset_name,folder,list_compound_id,prefix_name,trim)

def parse_args(args):
    parser = argparse.ArgumentParser(description='Collect .txt feature files into single .npy')

    parser.add_argument('--dataset_name', default='dataset', type=str)
    parser.add_argument('--compound_list_file', default='list_train.txt', type=str,                         help='Compound ids of a given dataset')
    parser.add_argument('--trim', action='store_true', default=False)

    args = parser.parse_args()
    return args
    
def cli_main():
    args = parse_args(sys.argv[1:])
    print(args)
    main(args)

if __name__ == "__main__":
    cli_main()
    print('End!')
