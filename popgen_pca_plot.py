#!/usr/bin/python

'''

popgen_pca_plot.py

Plotting pca following plink

WIP

'''
from sys import exit
from os import listdir, walk, symlink
from os.path import join, isfile, isdir

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def main(IN_EIGENVECS, IN_EIGENVALS, OUT_DIR):

    with open(IN_EIGENVECS) as input_eigenvecs, open(IN_EIGENVALS) as input_eigenvals:

        eigenvals = pd.read_csv(input_eigenvals, header = None)
        # print(eigenvals.describe)
        print('\n!!! EDIT TOTAL VAR FOR ACCURATE PLOT !!!\n')
        prop_var = pd.DataFrame()
        prop_var['val'] = eigenvals[0]
        prop_var['cumSum'] = eigenvals.cumsum()
        prop_var['prop'] = (prop_var['val'] / 42.65477 ) * 100
        prop_var['prop_cumSum'] = prop_var['prop'].cumsum()

        prop_of_var = prop_var['prop'].tolist()

        f, ax = plt.subplots(figsize=(8.5, 8.5))
        plot =sns.lineplot(data=eigenvals, x = eigenvals.index + 1, y=0, ax = ax)
        sns.scatterplot(data=eigenvals, x = eigenvals.index + 1, y=0, ax = ax)
        plt.xticks(eigenvals.index[::1] + 1)
        plt.xlabel("Principal Component")
        plt.ylabel("Eigenvalue")
        plot.figure.savefig(join(OUT_DIR, IN_EIGENVECS + '.scree.pdf'))

        # plot =sns.lineplot(data=prop_var, x = prop_var.index + 1, y='prop_cumSum', ax = ax)
        # sns.scatterplot(data=prop_var, x = prop_var.index + 1, y='prop_cumSum', ax = ax)
        # plt.xticks(eigenvals.index[::1] + 1)
        # plt.xlabel("Number of Principal Components")
        # plt.ylabel("Proportion of Variance explained")
        # plot.figure.savefig(join(OUT_DIR, IN_EIGENVECS + '.prop.pdf'))

        data_plot = pd.read_csv(input_eigenvecs, sep='\t')
        data_plot = data_plot.drop('#FID', axis=1)
        
        # data_plot['Location'] = data_plot['IID'].str.replace('_L001', '')
        # # data_plot['Location2'] = data_plot['Location'].str.replace('R0**', '')

        # locs_to_modify = data_plot['Location'].tolist()

        # new_locs = []

        # for samp in locs_to_modify:
        #     if len(samp.split('_')) > 1:
        #         new_locs.append(samp.split('_')[1])
        #     else:
        #         new_locs.append(samp)


        # ###All samples
        # refined_locs = []

        # with open(join(OUT_DIR, 'popfile.txt'), 'w') as out_file:
        #     for samp in new_locs:
        #         if samp.startswith('B'):
        #             if samp[7:].startswith('0'):
        #                 refined_locs.append(samp[8:])
        #                 out_file.write(samp[8:] + '\n')
        #             else:
        #                 refined_locs.append(samp[7:])
        #                 out_file.write(samp[7:] + '\n')
        #         else:
        #             refined_locs.append(samp)
        #             out_file.write(samp + '\n')


        # simplified_locs = []

        # with open(join(OUT_DIR, 'simple_popfile.txt'), 'w') as out_file:

        #     for samp in refined_locs:
        #         if 'SWDN' in samp:
        #             simplified_locs.append('SWDN')
        #             out_file.write('SWDN\n')
        #         elif 'UK' in samp.upper():
        #             simplified_locs.append('UK')
        #             out_file.write('UK\n')
        #         elif 'FRA' in samp:
        #             simplified_locs.append('FRA')
        #             out_file.write('FRA\n')
        #         elif samp == 'Ssx6':
        #             simplified_locs.append('UK')
        #             out_file.write('UK\n')

        #         else:
        #             print(samp)
        #             exit('EEK')

        ###SotonUK
        # simplified_locs = []
        # for samp in new_locs:
        #     # print (samp)
        #     if 'CpH' in samp:
        #         simplified_locs.append('Old')
        #         # out_file.write('SWDN\n')
        #     else:
        #         simplified_locs.append('New')

        ###LehFRA

        # simplified_locs = []
        # for samp in new_locs:
        #     print (samp)
        #     if 'CpC' in samp:
        #         simplified_locs.append('2013')
        #         # out_file.write('SWDN\n')
        #     else:
        #         simplified_locs.append('2014')

        individuals_with_Novo = ['VEM254', 'VEM252', 'VEM247', 'VEM239', 'VEM156', 'VEM154', 
                                 'VEM151', 'VEM140', 'VEM137', 'VEM098', 'VEM097', 'VEM089', 
                                 'VEM074', 'VEM068', 'VEM055', 'VEM053', 'VEM052', 'VEM051', 
                                 'VEM050', 'VEM044', 'VEM042', 'Blapi110', 'Blapi133', 'Blapi160',
                                 'Blapi166', 'Blapi261', 'Blapi29', 'Blapi327', 'Blapi356', 
                                 'Blapi377', 'Blapi412', 'Blapi47', 'Blapi497', 'Blapi531', 
                                 'Blapi540', 'Blapi557', 'Blapi628', 'Blapi96', 'EI002'] 
        
        # ### Creating plotting labels 
        data_plot['Location'] = np.where(data_plot['IID'].str.contains('VEM', regex=False), 'HISTORICAL', 'CONTEMPORARY')
        data_plot['Sequencer'] = np.where(data_plot['IID'].isin(individuals_with_Novo), 'NovaSeq', 'NextSeq')
        
        
        # data_plot = data_plot.drop(index=29)
        # data_plot = data_plot.drop(index=48)

        # print(data_plot.to_string())


        # first_three = data_plot.iloc[:, :4]
        first_five = data_plot.iloc[:, :6]
        first_five['Location'] = data_plot['Location']
        first_five['Sequencer'] = data_plot['Sequencer']
        # first_five['Location'] = data_plot['Location2']

        # first_three['Location'] = simplified_locs
        # print(first_three.describe)

        # first_five['Location'] = simplified_locs
        print(first_five.head())
        
        fig, ax = plt.subplots(figsize=(10, 10))

        # plot = sns.pairplot(first_three, hue='Location', diag_kind="hist", corner=True)
        # plot = sns.pairplot(first_five, hue='Location', diag_kind="None", diag_kws={"linewidth": 0, "shade": False}, corner=True)

        plot = sns.relplot(data=first_five, x='PC1',  y='PC2', hue='Location', style="Sequencer", legend=False) 
        # plot = sns.pairplot(first_three, diag_kind="hist", corner=True)

        #Label points
        for INDIV in first_five.IID:
            plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],first_five.IID[first_five.IID ==INDIV].to_string(index=False))


        #Label scots bees
        # for INDIV in ['VEM049', 'VEM050', 'VEM051', 'VEM052', 'VEM053', 'VEM055', 'VEM056']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],first_five.IID[first_five.IID ==INDIV].to_string(index=False))

        # for INDIV in ['VEM049', 'VEM050', 'VEM051', 'VEM053', 'VEM055', 'VEM056']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Scotland', fontsize=8, color='lightgrey')

        # for INDIV in ['VEM251']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Devon 1907', fontsize=8, color='lightgrey')
        
        # for INDIV in ['Blapi133', 'Blapi160']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Devon 2014', fontsize=8, color='lightgrey')

        # for INDIV in ['Blapi356']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Norfolk 2014', fontsize=8, color='lightgrey')

        # for INDIV in ['Blapi497']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Yorkshire 2014', fontsize=8, color='lightgrey')

        # for INDIV in ['Blapi110']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Kent 2014', fontsize=8, color='lightgrey')

        # for INDIV in ['VEM042', 'VEM098', 'VEM245']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Kent Museum', fontsize=8, color='lightgrey')

        # for INDIV in ['VEM239']:
        #     plt.text(first_five.PC1[first_five.IID==INDIV],first_five.PC2[first_five.IID==INDIV],'Gloustershire 1897', fontsize=8, color='lightgrey')


        # ### Hacky way of adding in the proportion of variance explained to the axis labels

        # replacements = {'PC1': 'PC1 (' + str(round(prop_of_var[0], 2)) + '%)',
        #                 'PC2': 'PC2 (' + str(round(prop_of_var[1], 2)) + '%)',
        #                 'PC3': 'PC3 (' + str(round(prop_of_var[2], 2)) + '%)',
        #                 'PC4': 'PC4 (' + str(round(prop_of_var[3], 2)) + '%)',
        #                 'PC5': 'PC5 (' + str(round(prop_of_var[5], 2)) + '%)'
        #                 }


        # for i in range(5):
        #     for j in range(5):
        #         if ax:
        #             try:
        #                 xlabel = plot.axes[i][j].get_xlabel()
        #                 ylabel = plot.axes[i][j].get_ylabel()
        #                 if xlabel in replacements.keys():
        #                     plot.axes[i][j].set_xlabel(replacements[xlabel])
        #                 if ylabel in replacements.keys():
        #                     plot.axes[i][j].set_ylabel(replacements[ylabel])
        #             except AttributeError:
        #                 pass
        # ###
        plt.legend().set_visible(False)
        # plot.savefig(join(OUT_DIR, IN_EIGENVECS + '.3.pairs.pdf'))
        # plot.savefig(join(OUT_DIR, '04_pc1_pc2.pdf'))
        plot.savefig(join(OUT_DIR,IN_EIGENVECS+'.pc1_pc2.edit.pdf'))

        print('Plot plotted')

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('plink_eigenvec', type=str, help='path to plink2 .eigenvec output file')
parser.add_argument('plink_eigenval', type=str, help='path to plink2 .eigenval output file')
# parser.add_argument('--list', dest='sp_list', type=str, help='if entered, prints plot for a subset defined by a list found at this path', default='NA')
parser.add_argument('-o', '--out_dir', dest='out_dir', type=str, help='path where figure should be written', default='./')
# parser.add_argument('out_dir', type=str, help='path to directory where the output should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.plink_eigenvec, args.plink_eigenval, args.out_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2021, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
