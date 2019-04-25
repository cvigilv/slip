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
        plt.hist(df[df[column] == state][Opt_SimilarityMeasure], density = True, color = colours[i], label = labels[i], alpha = 0.5, bins = bins)
    plt.xlabel(Opt_SimilarityMeasure)
    plt.ylabel('Relative count (%)')
    plt.legend(title = column)
    plt.title('{sim} distribution separated by {column}'.format(sim = Opt_SimilarityMeasure, column = column))
    plt.savefig(filename, dpi = 300)
    plt.cla()

# MAIN
print('SchuellerLab Ligand Priorization Pipepline - version {v}'.format(v = __version__))
print('Start time: {time}'.format(time = datetime.datetime.now()))

# Load output file to pandas dataframe
print('\nLoading output file to dataframe...')
in_df = pd.read_csv(Input_Interactions,
                    sep = '\t',
                    names = ['Fold', 'Query ligand ChEMBL ID', 'Hit target ChEMBL ID', Opt_SimilarityMeasure, 'Hit ligand ID', 'Query target ID', 'TP'],
                    header = None,
                    index_col = False) # Hardcoded options

print('--> Cleaning output file from buggy entries...',end = '\r')
in_df = in_df[in_df[Opt_SimilarityMeasure] >= 0]
print('--> Cleaning output file from buggy entries... DONE!')

if Output_Plots == True:
    print('--> Saving distribution plot for output file...',end = '\r')
    sns.boxplot(x = 'TP', y = Opt_SimilarityMeasure, data = in_df, palette = colours_TP)
    plt.ylabel(Opt_SimilarityMeasure)
    plt.title('{sim} distribution separated by TP'.format(sim = Opt_SimilarityMeasure))
    plt.savefig(Input_Interactions+'.distribution_boxplot.png', dpi = 300)
    plt.cla()

    TP_Histogram(in_df, 'TP', [0, 1], ['0', '1'], [color_neg, color_pos], Input_Interactions+'.distribution_hist.png')

    print('--> Saving distribution plot for output file... DONE!')
print('Loading output file to dataframe... DONE!')

# Filter output file based in similarity measure threshold
print('Filtering output by {sim} in range ]{min}, {max}]...'.format(sim = Opt_SimilarityMeasure, min = str(Opt_minSimilarity), max = str(Opt_maxSimilarity)))
in_df = in_df[(in_df[Opt_SimilarityMeasure] <= Opt_maxSimilarity) & (in_df[Opt_SimilarityMeasure] > Opt_minSimilarity)]

if Output_Plots == True:
    print('--> Saving distribution plot for filtered output file...',end = '\r')
    sns.boxplot(x = 'TP', y = Opt_SimilarityMeasure, data = in_df, palette = colours_TP)
    plt.ylabel(Opt_SimilarityMeasure)
    plt.title('{sim} distribution for range ]{min}, {max}] separated by TP'.format(sim = Opt_SimilarityMeasure,  min = str(Opt_minSimilarity), max = str(Opt_maxSimilarity)))
    plt.savefig(Input_Interactions+'.filtered_distribution_1.png', dpi = 300)
    plt.cla()

    TP_Histogram(in_df, 'TP', [0, 1], ['0', '1'], [color_neg, color_pos], Input_Interactions+'.filtered_dist_hist.png', bins = [Opt_minSimilarity + (Opt_maxSimilarity - Opt_minSimilarity)/20 * i for i in range(0,21)])

    print('--> Saving distribution plot for filtered output file... DONE!')
print('Filtering output by {sim} in range ]{min}, {max}]... DONE!'.format(sim = Opt_SimilarityMeasure, min = str(Opt_minSimilarity), max = str(Opt_maxSimilarity)))

# Add SMILES from SMILES file to query ligand based on ChEMBL ID
print('Adding SMILES for query ligand...', end = '\r')
smiles_df = pd.read_csv('/Users/cvigilv/Dropbox/Chembl22_goldStd3_max.txt.ul.co',
                        sep = '\t',
                        names = ['SMILES', 'Query ligand ChEMBL ID'],
                        header = None,
                        index_col = False)
in_df = pd.merge(in_df, smiles_df, on = 'Query ligand ChEMBL ID', how = 'left')
print('Adding SMILES for query ligand... DONE!' )

# Add information for query ligand and predicted target for next filters
print('Adding useful information to entries (Pfam of known targets, max clinical phase, number of atoms, etc.)')

broad_dict	= {}
pfam        = {}
max_phase   = {}
natoms      = {}
lig_pfam    = defaultdict(set)
targets     = defaultdict(set)

with open('/home/cvigilv/SLiP/Dependencies/chembl22_broad2.txt', 'r') as ints:
    for line in ints:
        tokens  = line.rstrip().split('\t')
        target  = tokens[0]
        ligid   = tokens[1]
        pfam_id = tokens[4]
        mphase  = tokens[11]

        natoms[ligid]       = tokens[12]
        pfam[target]        = pfam_id
        max_phase[ligid]    = mphase

        targets[ligid].add(target)
        lig_pfam[ligid].update(pfam_id.strip().split(','))      # Add PFam IDs of this target to the set of known IDs of this ligand


print('\t\t-> Loaded "chembl_broad.txt" to memory')



# print('Adding information for query ligand...')
# print('--> Loading broad selection of ChEMBL to memory...', end = '\r')
# broad_df = pd.read_csv(Input_Broad,
#                         sep = '\t',
#                         index_col = False,
#                         low_memory = False)
# broad_df = broad_df.rename(columns = {'chembl_id':'Target ChEMBL ID','chembl_id.1':'Ligand ChEMBL ID'})
# print('--> Loading broad selection of ChEMBL to memory... DONE!')
