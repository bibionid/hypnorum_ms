#!/usr/bin/python

'''

admixture_to_distruct.py

generate popfile for admixture - parse full sample names and create 'Early UK'

run with for i in {2..9}; do python admixture_to_distruct.py ${INPUT_VCF_NAME}.${i}.Q samples.txt ${i} ./ ; done

'''
# ##Load packages 
from os.path import join, isfile, isdir

import numpy  as np
import pandas as pd

import argparse


def main(IN_Q, IN_POP, K_IN, OUT_DIR):

    with open(IN_Q) as input_Q, open(IN_POP) as input_pop_labels:

        # ##Load input Q values from ADMIXTURE
        q_vals = pd.read_csv(input_Q, header = None, sep = ' ').astype(float)
        
        # ## Load the labels for the individuals
        pop_labs = pd.read_csv(input_pop_labels, header = None, names=['labs'])

        # ##Join the two things above 
        total = pd.concat([q_vals, pop_labs.reindex(q_vals.index)], axis=1)

        # ##Create the subset of Ealry UK samples
        early_UK = ['BOXCpH1PotdoUK', 'BOXCpH2HphUK', 'BOXCpH3SotonUK', 'BOXCpH4LdnUK', 'BOXCpH6SotonUK', 'BOXCpH8SotonUK','BOXCpH5HphUK', 'BOXCpH7SotonUK']

        # ##Making corrections to errors in labelling
        total['new_labs'] = np.where(total['labs'].str.len() > 28 , total['labs'].str.slice(start=19), total['labs'].str.slice(start=7)) #Trimming
        total['new_labs1'] = np.where(total['labs'].str.contains('SWDN', regex=False), 'SWDN', total['new_labs']) #Trimming
        total['new_labs2'] = np.where(total['new_labs1'].str.contains('Ssx6', regex=False), 'SsxUK', total['new_labs1']) # Correcting Ssx6 mislabel
        total['new_labs3'] = np.where(total['new_labs2'].str.contains('NcleUk', regex=False), 'NcleUK', total['new_labs2']) # Correcting NcleUk mislabel
        total['new_labs4'] = np.where(total['new_labs3'].str.contains('leFRA', regex=False), 'LilleFRA', total['new_labs3']) # Correcting NcleUk mislabel
        total['new_labs5'] = np.where(total['new_labs4'].str.startswith('0'), total['new_labs4'].str.slice(start=1), total['new_labs4']) # Remove odd 0
        total['new_labs6'] = np.where(total['labs'].isin(early_UK), 'EarlyUK', total['new_labs5']) # Relable the Early UK samples

        #Create an order for the samples based on geographic location
        geog_order = {'SWDN': 0, 'BdxFRA': 1, 'AgrsFRA': 2, 'LehFRA': 3, 'LilleFRA': 4,
            'EarlyUK': 5, 'SotonUK': 6, 'PlyUK': 7, 'SsxUK': 8, 'CardUK': 9, 'CambsUK': 10,
            'BgrUK': 11, 'McrUK': 12, 'HullUK': 13, 'NcleUK': 14}

        # ##Create the reordering variable and sort by it
        total['geog_rank'] = total['new_labs6'].map(geog_order)
        total = total.sort_values(by = ['geog_rank'], ascending=False)

        # ##Uncomment to see the output of this process
        # print(total.to_string())
        # print()

        #Create the input files for distruct2.3
        with open(join(OUT_DIR, 'popfile_EarlyUK_'+ str(K_IN)), 'w') as out_file:
            total['new_labs6'].to_csv(out_file, header=None, index=None, sep=' ', mode='a')

        with open(join(OUT_DIR, IN_Q.replace('.Q', '') + '_' + str(K_IN)) + '.Q.geogSorted.EarlyUK', 'w') as out_file:
            total[range(K_IN)].to_csv(out_file, header=None, index=None, sep=' ', mode='w', float_format='%.6f')


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('Q_in', type=str, help='path to .Q file from admixture')
parser.add_argument('popfile_in', type=str, help='path to popfile generated by pca script')
parser.add_argument('k_in', type=int, help='number of k from admixture analysis')
parser.add_argument('out_dir', type=str, help='path where output files should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.Q_in, args.popfile_in, args.k_in, args.out_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2024, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
