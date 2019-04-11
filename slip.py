#!/usr/bin/env python
# -*- coding: utf-8 -*

__version__ = 0.0

# DEPENDENCIES
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from configparser import ConfigParser
import datetime


# SETTINGS
configs = ConfigParser()
configs.read(sys.argv[1])

Input_Interactions = configs.get('Input', 'Interactions file')
Input_Broad        = configs.get('Input', 'Broad file')
Input_ChEMBL       = configs.get('Input', 'ChEMBL version')

Opt_JustTarget       = configs.get('Options', 'Keep target')
Opt_JustLigand       = configs.get('Options', 'Keep ligand')
Opt_PfamCutoff       = configs.getint('Options', 'Pfam cutoff')
Opt_SimilarityMeasure= configs.get('Options', 'Similarity measure used')
Opt_TopX             = configs.getint('Options', 'Top X entries')
Opt_minSimilarity    = configs.getfloat('Options', 'min(Similarity)')
Opt_maxSimilarity    = configs.getfloat('Options', 'max(Similarity)')
Opt_minClinicalPhase = configs.getint('Options', 'min(Clinical Phase)')
Opt_maxClinicalPhase = configs.getint('Options', 'max(Clinical Phase)')

Output_Directory = configs.get('Output', 'Output directory')
Output_File      = configs.get('Output', 'Output file')
Output_Plots     = configs.getboolean('Output', 'Generate plots')

if Opt_SimilarityMeasure == '': Opt_SimilarityMeasure = 'Similarity measure'

color_neg = '#527AB2'
color_pos = '#FF4528'
colours_TP = {0: color_neg, 1: color_pos}


# FUNCTIONS
def TP_Histogram(df, column, states, labels, colours, filename, bins = [0.05*i for i in range(0,21)]):
    for i, state in enumerate(states):
        plt.hist(df[df[column] == state]['Similarity measure'], density = True, color = colours[i], label = labels[i], alpha = 0.5, bins = bins)
    plt.xlabel(Opt_SimilarityMeasure)
    plt.ylabel('Relative count (%)')
    plt.legend(title = column)
    plt.title('{sim} distribution separated by {column}'.format(sim = Opt_SimilarityMeasure, column = column))
    plt.savefig(filename, dpi = 300)
    plt.cla()

# MAIN
print('SchuellerLab Ligand Priorization Pipepline - version {v}'.format(v = __version__))
print('Start time: {time}'.format(time = datetime.datetime.now()))
print('\nLoading output file to dataframe...',end = '\r')
in_df = pd.read_csv(Input_Interactions,
                    sep = '\t',
                    names = ['Fold', 'Query ligand ChEMBL ID', 'Hit target ChEMBL ID', Opt_SimilarityMeasure, 'Hit ligand ID', 'Query target ID', 'TP'],
                    header = None,
                    index_col = False) # Hardcoded options

print('Loading output file to dataframe... DONE!')
print('\nSummary of output file:')
print(in_df.info(memory_usage='deep'))

if Output_Plots == True:
    print('Saving distribution plot for output file...',end = '\r')
    sns.boxplot(x = 'TP', y = 'Similarity measure', data = in_df, palette = colours_TP)
    plt.ylabel(Opt_SimilarityMeasure)
    plt.title('{sim} distribution separated by TP'.format(sim = Opt_SimilarityMeasure))
    plt.savefig(Input_Interactions+'.distribution_boxplot.png', dpi = 300)
    plt.cla()

    TP_Histogram(in_df, 'TP', [0, 1], ['0', '1'], [color_neg, color_pos], Input_Interactions+'.distribution_hist.png')

    print('Saving distribution plot for output file... DONE!')


print('\nFiltering output by {sim} in range ]{min}, {max}]...'.format(sim = Opt_SimilarityMeasure, min = str(Opt_minSimilarity), max = str(Opt_maxSimilarity)), end = '\r')
in_df = in_df[(in_df['Similarity measure'] <= Opt_maxSimilarity) & (in_df['Similarity measure'] > Opt_minSimilarity)]
print('Filtering output by {sim} in range ]{min}, {max}]... DONE!'.format(sim = Opt_SimilarityMeasure, min = str(Opt_minSimilarity), max = str(Opt_maxSimilarity)))
print(in_df.info(memory_usage = 'deep'))

if Output_Plots == True:
    print('Saving distribution plot for filtered output file...',end = '\r')
    sns.boxplot(x = 'TP', y = 'Similarity measure', data = in_df, palette = colours_TP)
    plt.ylabel(Opt_SimilarityMeasure)
    plt.title('{sim} distribution for range ]{min}, {max}] separated by TP'.format(sim = Opt_SimilarityMeasure,  min = str(Opt_minSimilarity), max = str(Opt_maxSimilarity)))
    plt.savefig(Input_Interactions+'.filtered_distribution_1.png', dpi = 300)
    plt.cla()

    TP_Histogram(in_df, 'TP', [0, 1], ['0', '1'], [color_neg, color_pos], Input_Interactions+'.filtered_dist_hist.png', bins = [Opt_minSimilarity + (Opt_maxSimilarity - Opt_minSimilarity)/20 * i for i in range(0,21)])

    print('Saving distribution plot for filtered output file... DONE!')
