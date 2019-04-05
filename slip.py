#!/usr/bin/env python
# -*- coding: utf-8 -*

__version__ = 0.0

# DEPENDENCIES
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# SETTINGS
config = ConfigParser()
config.read(sys.argv[2])

Input_Interactions = configs.get('Input', 'Interactions file')
Input_Broad        = configs.get('Input', 'Broad file')
Input_ChEMBL       = configs.get('Input', 'ChEMBL version')

Opt_JustTarget       = configs.get('Options', 'Keep target')
Opt_JustLigand       = configs.get('Options', 'Keep ligand')
Opt_PfamCutoff       = configs.getint('Options', 'Pfam cutoff')
Opt_TopX             = configs.getint('Options', 'Top X entries')
Opt_minSimilarity    = configs.getint('Options', 'min(Similarity)')
Opt_maxSimilarity    = configs.getint('Options', 'max(Similarity)')
Opt_minClinicalPhase = configs.getint('Options', 'min(Clinical Phase)')
Opt_maxClinicalPhase = configs.getint('Options', 'max(Clinical Phase)')

Output_Directory = configs.get('Output', 'Output directory')
Output_File      = configs.get('Output', 'Output file')
Output_Plots     = configs.get('Output', 'Generate plots')

# FUNCTIONS

# MAIN
print('SchuellerLab Ligand Priorization Pipepline - version {v}'.format(v = __version__))


in_df = pd.read_csv(Input_Interactions,
                    sep = '\t',
                    names = ['Fold', 'Query ligand ChEMBL ID', 'Hit target ChEMBL ID', 'Similarity measure', 'Hit ligand ID', 'Query target ID', 'TP'],
                    header = False,
                    index_col = False) # Hardcoded options
