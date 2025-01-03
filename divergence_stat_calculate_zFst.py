#!/usr/bin/python

'''

divergence_stat_calculate_zFst.py

calculate window-wise zFst using the output of pixy

window zFst = (window Fst - mean Fst across all windows) / std dev of mean Fst across all windows

'''

from os.path import join
from itertools import permutations

import pandas as pd
import argparse

def main(PIXY_IN, OUT_DIR):

    #Load pixy Fst ouput 
    pixyFst_df = pd.read_csv(PIXY_IN, delimiter='\t', header=0, names=['pop1','pop2','chromosome','window_pos_1','window_pos_2','avg_wc_fst','no_snps'])

    #filter windows with < 10 SNPs
    pixyFst_df = pixyFst_df[pixyFst_df['no_snps'] >= 10]

    #calculate mean Fst and sd across all windows
    pixyFst_df['meanFst'] = pixyFst_df['avg_wc_fst'].mean()
    pixyFst_df['meanFst_sd'] = pixyFst_df['avg_wc_fst'].std()

    #calculate zFst
    pixyFst_df['zFst'] = ( pixyFst_df['avg_wc_fst'] - pixyFst_df['meanFst'] ) / pixyFst_df['meanFst_sd']

    #save output     
    pixyFst_df[['chromosome', 'window_pos_1', 'window_pos_2', 'zFst']].to_csv(join(OUT_DIR, PIXY_IN.replace('.txt', '.zFst.txt')), sep='\t', index=False, header=False)
    

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('pixy_in', type=str, help='path to dir containing per scaffold pixy results')
parser.add_argument('out_dir', type=str, help='path where output files should be written')
# parser.add_argument('-s', '--scaffold', type=str, help='option to only focus on one scaffold and plot a line plot', nargs='?', const='NA', default='NA')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.pixy_in, args.out_dir)


__author__ = "Will Nash"
__copyright__ = "Copyright 2024, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
