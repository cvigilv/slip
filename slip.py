'''
SLIP - SchuellerLab LIgand Pipeline
version 0.2

Authors:	Carlos Vigil, Andreas Schüller

To run:
    > python3 slip.py "input_file" "configuration_file"
    where       input_file is a .txt file, configuration_file is a .py file;
                both have specific formatting.

Changelog:
    -   2018.06.08	0.1     Carlos	First version
    -	2018.11.14	0.2		Carlos	Added alternative directory option for missing SMILES


Dependencies:
    -   Python libraries:
            ·   pymysql
            ·   prettytable
    -   MySQL:
            ·   chembl_23
    -   MOE 2018
'''
#!/usr/bin/env python
# -*- coding: utf-8 -*


# LIBRARY
import sys
import datetime
import os
import pymysql
import prettytable
import subprocess
from collections    import defaultdict
from operator       import itemgetter
from configparser   import ConfigParser

# CONFIG LOADING
config = ConfigParser()

config.read(sys.argv[2])

Keep_Predicted_Target   = config.get('options','Keep_Predicted_Target')
commonPfam_cutoff       = config.getint('options','commonPfam_cutoff')
TopX_cutoff             = config.getint('options','TopX_cutoff')
Tc_min                  = config.getfloat('options','Tc_min')
Tc_max                  = config.getfloat('options','Tc_max')
mphase_min              = config.getint('options','mphase_min')
mphase_max              = config.getint('options','mphase_max')
ALTsmiles		        = config.get('options','alt_smiles_path')
ALTbroad				= config.get('options','alt_broad_path')


# VARIABLES
list_fileName   = []
list_lineCount  = []
list_filterID   = []

#   FUNCTIONS
#+  OTHER FUNCTIONS
def progressPercentage(currentEntry):
    '''
    Calculate the percentage of progress from an specified step
    '''
    global list_lineCount

    total   = int(list_lineCount[-1])
    Percent = round( ((currentEntry/total)*100), 1)

    return str(Percent)+' / 100%'

def initializeSLIP():
    '''
    Initialize SLiP Pipeline:
    i.   Creates directory for output files
    ii.  Prints out file to process, output folder name, amount of entries in input file and SLiP configuration.
    '''

    global cvout_file, config_file, fileName, list_fileName, list_lineCount, Keep_Predicted_Target

    fileName         = cvout_file.split('/')[-1]; list_fileName.append(fileName)
    folderExtension  = ('-'+datetime.datetime.now().strftime("%Y.%m.%d_%H:%M:%S")).rstrip()
    folderName       = fileName+folderExtension
    configString     = 'Pfam cutoff == {} - Top {} - Tc == [{},{}[ - MaxPhase == [{},{}['.format(commonPfam_cutoff, TopX_cutoff,Tc_min,Tc_max,mphase_min,mphase_max)
    if Keep_Predicted_Target != '': 'Keep only entries with "{}" as Predicted Target - {}'.format(Keep_Predicted_Target,configString)

    print(">> Input information <<\n")
    print("input_file   ->\t\t{}".format(os.path.abspath(cvout_file)))
    print("out_folder   ->\t\t{}".format(os.path.abspath(folderName)))

    os.system('mkdir {1}; mv {0} {1}/; mv {2} {1}/'.format(cvout_file, folderName, config_file))
    os.chdir('{}/'.format(folderName))

    lines_cvout = os.popen('wc -l {}'.format(cvout_file)).read().split(' ')[0]
    list_lineCount.append(lines_cvout)

    print("input_lines  ->\t\t{:,}".format(int(lines_cvout)))
    print("\nconfig_file  ->\t\t{}".format(os.path.abspath(config_file)))
    print("configs      ->\t\t{}".format(configString))

    print("\nChanging work directory to {}".format(os.getcwd()))

#+  PIPELINE FUNCTIONS
def Clean_Input(cvout_file):
    '''
    Cleans input file to eliminate posible errors:
    i.   Ignore cases with Tc == -99 (False Negatives)
    ii.  Ignore cases with TP == 1.0 (True Positives)
    '''

    global fileName, list_fileName, list_lineCount

    print("\n i.\tCleaning input file")

    fileName            = fileName+'.clean'; list_fileName.append(fileName)
    outFile             = open(fileName,'a+')
    entryN              = 0
    dismissedEntries    = 0

    with open(cvout_file,'r') as file:
        for line in file:
            tokens = line.rstrip().split('\t')

            Tc = float(tokens[3])
            TP = float(tokens[6])

            entryN += 1
            print ("\t\t-> {}".format(progressPercentage(entryN)), end='\r')

            # Check "Tc" & "TP", if not problematic then write line to output file
            if Tc != -99:
                if TP == 0.0:
                    outFile.write(line)
                else:
                    dismissedEntries +=1
            else:
                dismissedEntries +=1

    outFile.close()

    lines_clean = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
    list_lineCount.append(lines_clean)

    print ("\t\t-> Amount of predictions kept:\t\t\t{:,}\n\t\t-> Amount of predictions dismissed:\t\t{:,}".format(int(lines_clean),dismissedEntries))

def Select_Target(cvout_file):
    '''
    Keep specified target in "slip_conf.py" file:
    i.  Extract variable from "slip_conf.py" file.
    ii. Check if "Predicted_Target" is the same as the specified target. If TRUE, then write line.
    #   Currently working only for CHEMBL targets.
    '''
    global fileName, list_fileName, list_lineCount

    if Keep_Predicted_Target != '':
        if 'CHEMBL' in Keep_Predicted_Target:

            keepTarget   = Keep_Predicted_Target
            fileName     = fileName+'.'+keepTarget; list_fileName.append(fileName)

            print(' ii.\tKeeping entries for predicted target "{}"'.format(Keep_Predicted_Target))
            os.system("awk '($3 == \"{0}\")' {1} >> {2}".format(keepTarget,list_fileName[-2],fileName))

            lines_target = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
            list_lineCount.append(lines_target)

            print ("\t\t-> Amount of predictions for specified target:\t{:,}".format(int(lines_target)))

        else:
            print("[FATAL ERROR : Invalid specified predicted target ID, check 'slip_conf.py' file]")
            exit(1)

def Add_SMILES(cvout_file):
    '''
    Add canonical SMILES to each entry:
    i.  Load to memory all the canonical SMILES from file "Chembl22_goldStd3_max.txt.ul.co" and of alternative SMILES file
    ii.
    '''
    global fileName, list_fileName, list_lineCount

    print(' iii.\tAdding canonical SMILES')

    fileName = fileName+'.smiles'; list_fileName.append(fileName)
    out      = open(fileName,'a+')

    smiles   	= {}
    co			= {}
    alt_smiles 	= {}
    alt_co		= {}
    entry    	= 0

    with open('/home/cvigilv/SLiP/Dependencies/Chembl22_goldStd3_max.txt.ul.co') as file:
        for line in file:
            tokens              = line.rstrip().split('\t')
            smiles[tokens[1]]   = tokens[0]
            co[tokens[2]]		= tokens[1]

    if ALTsmiles not in [None, 'None']:
        with open(ALTsmiles) as file:
            for line in file:
                tokens              = line.rstrip().split('\t')
                alt_smiles[tokens[1]]   = tokens[0]
                alt_co[tokens[1]]		= tokens[2]

    with open(cvout_file) as file:
        for line in file:
            tokens              = line.rstrip().split('\t')
            chemblid            = tokens[1]

            entry += 1
            print ("\t\t-> {}".format(progressPercentage(entry)),end = '\r')

            if chemblid in smiles and tokens[5] in smiles:
                out.write(line.rstrip() + '\t' + smiles[chemblid] + '\t' + smiles[tokens[5]] + '\n')
            elif chemblid in alt_smiles and tokens[5] in smiles:
                out.write(line.rstrip() + '\t' + alt_smiles[chemblid] + '\t' + smiles[tokens[5]] + '\n')
            elif chemblid in smiles and tokens[5] in alt_smiles:
                out.write(line.rstrip() + '\t' + smiles[chemblid] + '\t' + alt_smiles[tokens[5]] + '\n')
            elif chemblid in alt_smiles and tokens[5] in alt_smiles:
                out.write(line.rstrip() + '\t' + alt_smiles[chemblid] + '\t' + alt_smiles[tokens[5]] + '\n')
            else:
                sys.stderr.write('ERROR: SMILES not found for ID %s\n' % chemblid)

    out.close()

    lines_smiles = os.popen('wc -l {}'.format(cvout_file+'.smiles')).read().split(' ')[0]
    list_lineCount.append(lines_smiles)

def Add_Info(cvout_file):
    '''
    Add relevant information for filtering steps:
    i.   Load to memory "Query_Ligand" Pfam, known targets, max clinical phase & number of atoms.
    ii.  Load to memory "Hit_Ligand" Pfam and known targets.
    iii.

    '''
    global fileName, list_fileName, list_lineCount

    print('\n iv.\tAdding useful information to entries (Pfam of known targets, max clinical phase, number of atoms, etc.)')

    fileName    = fileName+'.pfam+mphase+natoms'; list_fileName.append(fileName)
    out         = open(fileName,'a+')
    entryN      = 0

    # Information loaded to memory is stored here
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

    with open(cvout_file, 'r') as cvout:
        for line in cvout:
            tokens      = line.rstrip().split('\t')
            qry_ligid   = tokens[1]                                 # Query ligand
            pred_target = tokens[2]                                 # Predicted target
            #         Original line          Known targets of the query ligand     PFam IDs of the known trgs   PFams of predicted trg Max clinical phase of query
            out.write(line.rstrip() + '\t '+ ','.join(targets[qry_ligid]) + '\t' + ','.join(lig_pfam[qry_ligid]) + '\t' + pfam[pred_target] + '\t' + max_phase[qry_ligid] + '\t' + natoms[qry_ligid] + '\n')#  ''+'\t'+ min_act[(qry_ligid, target)] + '\t' + avg_act[(qry_ligid, target)] + '\t' + max_act[(qry_ligid, target)] + '\n')

            entryN += 1
            print ("\t\t-> {}".format(progressPercentage(entryN)),end = '\r')
    out.close()

    lines_pfam = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
    list_lineCount.append(lines_pfam)

def MySQL(cvout_file):
    '''
    Check if prediction exist in chembl_23 database:
    i.  Get information from input file and assign to variables.
    ii. Run MySQL query to check existance of ligand-protein interaction. If FALSE, write line to output file.
    #   This segment has issues with some harcoded variables.
    #+  Slowest step in the pìpeline. Need to see way to minimize this issue.
    '''
    global fileName, list_fileName, list_lineCount, commonPfam_cutoff

    print('\n v.\tMySQL query filter')


    fileName = fileName+'.MySQL'; list_fileName.append(fileName)
    out      = open(fileName,'a+')
    entryN   = 0

    db = pymysql.connect("localhost","root","123","chembl_23")

    with open(cvout_file, 'r') as ints:
        for line in ints:
            tokens  = line.rstrip().split('\t')
            pfam    = tokens[10].strip().split(',')             # Pfam IDs of the predicted target
            pfam    = [x for x in pfam if x.startswith('PF')]

            if len(pfam) == 0:
                continue

            ligid   = tokens[2]
            target  = tokens[3]
            cursor  = db.cursor()
            results = []

            sql = """SELECT * FROM activities act
                JOIN molecule_dictionary    AS md ON act.molregno    = md.molregno
                JOIN assays                 AS a  ON a.assay_id      = act.assay_id
                JOIN target_dictionary      AS td ON a.tid           = td.tid
                LEFT JOIN target_components AS tc ON td.tid          = tc.tid
                LEFT JOIN component_domains AS cd ON tc.component_id = cd.component_id
                LEFT JOIN domains           AS do ON cd.domain_id    = do.domain_id
                WHERE md.chembl_id = '%s' AND (td.chembl_id = '%s' OR (do.domain_type = 'Pfam-A' AND do.source_domain_id IN ('%s')))""" % (ligid, target, "','".join(pfam))

            cursor.execute(sql)

            entryN += 1
            print ("\t\t-> {}".format(progressPercentage(entryN)),end = '\r')

            if cursor.rowcount <= 0:                            # Currently hardcoded, eventhough the configuration variable exist. CHECK FOR NEXT VERSION!
                out.write(line)

    db.close()
    out.close()

    lines_mysql = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
    list_lineCount.append(lines_mysql)

def Add_Molport(cvout_file,instock):
    '''
    Add MolPort ID for experimental checking of predictions:
    i.  Load to memory MolPort db ("InStock" or "Boutique").
    ii. Write all MolPort ID's (concatenated) for each "Query_Ligand" of input file.
    #   2 columns are added in the pipeline: "InStock" Molport & "Boutique" Molport.
    #+  Currently (v_1.0) only adding "InStock" Molport due too memory issues in archimedes.
    '''
    global fileName, list_fileName, list_lineCount

    print('\n vi.\tAdding MolPort ID\'s')

    entryN = 0
    chembl2molport = defaultdict(set)

    if instock == 1:
        fileName    = fileName+'.instock'; list_fileName.append(fileName)
        out         = open(fileName,'a')

        with open('/home/cvigilv/SLiP/Dependencies/molport.instock.can', 'r') as file:
            for line in file:
                line    = line.rstrip('\n')
                tokens  = line.split('\t')
                smiles  = tokens[0]
                molport = tokens[1]

                chembl2molport[smiles].add(molport)

    elif instock == 0:
        fileName    = fileName+'.boutique'; list_fileName.append(fileName)
        out         = open(fileName,'a')

        with open('/home/cvigilv/SLiP/Dependencies/molport.notinstock.can', 'r') as file:
            for line in file:
                line    = line.rstrip('\n')
                tokens  = line.split('\t')
                smiles  = tokens[0]
                molport = tokens[1]

                chembl2molport[smiles].add(molport)

    with open(cvout_file, 'r') as cvout:
        for line in cvout:
            line    = line.rstrip()
            tokens  = line.split('\t')
            smiles  = tokens[0]

            entryN += 1
            print ("\t\t-> {}".format(progressPercentage(entryN)),end = '\r')

            if smiles in chembl2molport:
                out.write(line + '\t' + ','.join(list(chembl2molport[smiles])) + '\n')
                out.flush()                     # Problems with size of line solved with this.
            else:
                out.write(line + '\t\n')
                out.flush()                     # Problems with size of line solved with this.

    out.close()

    lines_molport    = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
    list_lineCount.append(lines_molport)

def Tc_Cutoff(cvout_file):
    '''
    Keep entries with "Tc" greater or equal to "x" and less than "y":
    i. Checks if "Tc" of each line is between desired range specified in the configuration file.
    '''
    global fileName, list_fileName, list_lineCount, Tc_min, Tc_max

    print('\n vii.\tKeeping cases with Tc in [{},{}[ range'.format(Tc_min,Tc_max))

    Tc_name     = '.Tc_ge{}_lt{}'.format(Tc_min,Tc_max)
    fileName    = fileName+Tc_name; list_fileName.append(fileName)

    outFile = open(fileName,'a')
    Tc_min  = float(Tc_min)
    Tc_max  = float(Tc_max)
    entryN  = 0

    with open(cvout_file) as file:
        for line in file:
            tokens  = line.rstrip().split('\t')
            Tc      = float(tokens[4])

            entryN += 1
            print ("\t\t-> {}".format(progressPercentage(entryN)),end = '\r')

            if Tc < Tc_max:
                if Tc >= Tc_min:
                    outFile.write(line)


    outFile.close()

    lines_Tc = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
    list_lineCount.append(lines_Tc)

def MaxPhase_Cutoff(cvout_file):
    '''
    Keep entries with "Max_Phase" greater or equal to "x" and less than "y":
    i. Checks if "MaxPhase" of each line is between desired range specified in the configuration file.
    '''
    global fileName, list_fileName, list_lineCount, mphase_min, mphase_max

    print('\n viii.\tKeeping cases with "MaxPhase" in [{},{}[ range'.format(mphase_min,mphase_max))

    maxPhase_name   = '.maxPhase_ge{}_lt{}'.format(Tc_min,Tc_max)
    fileName        = fileName+maxPhase_name; list_fileName.append(fileName)

    outFile         = open(fileName,'a')
    maxPhase_min    = float(mphase_min)
    maxPhase_max    = float(mphase_max)
    entryN          = 0

    with open(cvout_file) as file:
        for line in file:
            tokens  = line.rstrip().split('\t')
            maxPhase= float(tokens[12])

            entryN += 1
            print ("\t\t-> {}".format(progressPercentage(entryN)),end = '\r')

            if maxPhase < maxPhase_max:
                if maxPhase >= maxPhase_min:
                    outFile.write(line)


    outFile.close()

    lines_MaxPhase = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
    list_lineCount.append(lines_MaxPhase)

def TopX(cvout_file):
    '''
    Keeps top "N" (specified in configuration file) predictions:
    i.  Load input to memory and enumerate entries.
    ii. Sort entries based in "Tc"
    iii.Write first "N" entries to output file.
    '''
    global fileName, list_fileName, list_lineCount, TopX_cutoff

    print('\n ix.\tKeeping Top {} cases'.format(TopX_cutoff))


    fileName    = fileName+'.top'+str(TopX_cutoff); list_fileName.append(fileName)
    outFile     = open(fileName,'a')
    N_cutoff    = int(TopX_cutoff)

    # Variable to load to memory
    cvout = []

    with open(cvout_file) as file:
        for line in file:
            tokens = line.rstrip().split('\t')
            cvout.append(tokens)

    # Sort input file based in "Tc"
    cvout.sort(key=lambda k: float(k[4]), reverse=True)

    for entryN,line in enumerate(cvout):
        if entryN == N_cutoff:
            break

        outFile.write('\t'.join(line)+'\n')
        print ("\t\t-> {}".format(progressPercentage(entryN)),end = '\r')


    outFile.close()

    lines_TopX = os.popen('wc -l {}'.format(fileName)).read().split(' ')[0]
    list_lineCount.append(lines_TopX)

# MAIN
#+ Assign command arguments to variables for easier lecture
cvout_file = sys.argv[1]
config_file = sys.argv[2]

#+ Initialize SLiP
os.system('clear')
print('\nSLIP - SchuellerLab LIgand Pipeline\n%s\n' % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
initializeSLIP()

#+ Filtering Stage
print("\n\n>> SLiP - Ligand Filtering <<")

Clean_Input(fileName)
Select_Target(fileName)
Add_SMILES(fileName)
Add_Info(fileName)
MySQL(fileName)
Add_Molport(fileName,1)
#Tc_Cutoff(fileName)
MaxPhase_Cutoff(fileName)
TopX(fileName)

print("-> Ligand Filtering : DONE!\n")

#+ Create table with results of "Ligand Filtering" section
log = open('slip.log','w+')

if (list_fileName and list_lineCount):
    Results = prettytable.PrettyTable(["File name", "Entries"])

    Results.align["File name"]  = "l"           #++ Align text to left side
    Results.align["Entries"]    = "r"           #++ Align text to right side

    a = 0
    for i in list_fileName:
        Results.add_row([i,'{:,}'.format(int(list_lineCount[a]))])
        a += 1

    print('\nTable 1: Amount of entries for each step of SLiP "Ligand Filtering"')
    print(Results)

    # log.write('SLIP - Log')
    # log.write('\nTable 1: Amount of entries for each step of SLiP "Ligand Filtering"')
    # log.write(str(Results))

#+ Preparation Stage
print("\n\n>> SLiP - Ligand Preparation<<\n")
os.system('moebatch -load /home/cvigilv/SLiP/Dependencies/slip_prepare.svl -exec "PrepareLigands \'{}\'"'.format(fileName))
