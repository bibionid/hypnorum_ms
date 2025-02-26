#!/usr/bin/python

'''

divergence_stat_calculate_zFst.py

calculate window-wise zFst using the output of pixy

window zFst = (window Fst - mean Fst across all windows) / std dev of mean Fst across all windows

'''

from os.path import join

import pandas as pd
import argparse
from numpy import sqrt

from sys import exit

def main(PIXY_IN, OUT_DIR):

    # ###Load pixy Fst ouput 
    pixyFst_df = pd.read_csv(PIXY_IN, delimiter='\t', header=0, names=['pop1','pop2','chromosome','window_pos_1','window_pos_2','avg_wc_fst','no_snps'])

    # ###add pop pair column
    pixyFst_df['pop_pair'] = pixyFst_df['pop1'] + '_' + pixyFst_df['pop2']

    # ###filter windows with < 10 SNPs
    pixyFst_df = pixyFst_df[pixyFst_df['no_snps'] >= 10]

    # ###calculate mean Fst and sd across all windows for each population pair
    pop_pair_means = pixyFst_df.groupby('pop_pair')['avg_wc_fst'].mean()
    pop_pair_std = pixyFst_df.groupby('pop_pair')['avg_wc_fst'].std()

    # ###calculate zFst based on per population pair specific mean and sd
    pixyFst_df['zFst'] = ( pixyFst_df['avg_wc_fst'] - pop_pair_means[pixyFst_df['pop_pair']].values ) / pop_pair_std[pixyFst_df['pop_pair']].values

    #save output     
    pixyFst_df[['pop_pair', 'chromosome', 'window_pos_1', 'window_pos_2', 'zFst']].to_csv(join(OUT_DIR, PIXY_IN.replace('.txt', '.pairwise.zFst.txt')), sep='\t', index=False, header=False)

    # ###Part 2 creating a per bin average zFst for inter and intra group comparisons

    # ###create in and out groups
    out_pairs = ['AgrsFRA_HullUK', 'AgrsFRA_McrUK', 'AgrsFRA_NcleUK', 'HullUK_BdxFRA', 'HullUK_LehFRA', 'HullUK_LilleFRA', 'McrUK_BdxFRA', 'McrUK_LehFRA', 'McrUK_LilleFRA', 'NcleUK_BdxFRA', 'NcleUK_LehFRA', 'NcleUK_LilleFRA']
    in_pairs  = ['AgrsFRA_BdxFRA', 'AgrsFRA_LehFRA', 'AgrsFRA_LilleFRA', 'BdxFRA_LehFRA', 'BdxFRA_LilleFRA', 'LehFRA_LilleFRA', 'McrUK_HullUK', 'NcleUK_HullUK', 'NcleUK_McrUK']

    # ###subset out the in and out groups and add a label column
    out_set = pixyFst_df[pixyFst_df['pop_pair'].isin(out_pairs)]
    out_set['label'] = 'out_group'
    in_set  = pixyFst_df[pixyFst_df['pop_pair'].isin(in_pairs)]
    in_set['label'] = 'in_group'
    
    # ###add a per window mean to each set
    out_set['mean_zFst'] = out_set.groupby(['chromosome', 'window_pos_1', 'window_pos_2'])['zFst'].transform('mean')
    out_set['mean_zFst_sd'] = out_set.groupby(['chromosome', 'window_pos_1', 'window_pos_2'])['zFst'].transform('std')
    out_set['mean_zFst_se'] = out_set['mean_zFst_sd']/ sqrt(out_set.groupby(['chromosome', 'window_pos_1', 'window_pos_2'])['zFst'].transform('count'))
    in_set['mean_zFst'] = in_set.groupby(['chromosome', 'window_pos_1', 'window_pos_2'])['zFst'].transform('mean')
    in_set['mean_zFst_sd'] = in_set.groupby(['chromosome', 'window_pos_1', 'window_pos_2'])['zFst'].transform('std')
    in_set['mean_zFst_se'] = in_set['mean_zFst_sd']/ sqrt(in_set.groupby(['chromosome', 'window_pos_1', 'window_pos_2'])['zFst'].transform('count'))

    # ###concatenate the two sets
    zFst_df = pd.concat([out_set, in_set])      

    # ###Create the second output file
    zFst_df[['pop_pair', 'chromosome', 'window_pos_1', 'window_pos_2', 'label', 'mean_zFst', 'mean_zFst_se']].to_csv(join(OUT_DIR, PIXY_IN.replace('.txt', '.windowMeans.zFst.txt')), sep='\t', index=False, header=False)

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('pixy_in', type=str, help='path to dir containing per scaffold pixy results')
parser.add_argument('out_dir', type=str, help='path where output files should be written')

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
