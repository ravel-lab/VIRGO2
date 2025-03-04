#######VIRGO2
#######Authors: Michael France, Issac Chaudry
#######Contact: mfrance@som.umaryland.edu

####DEPENDENCIES required to be in PATH
###bowtie2
###samtools

"""
Copyright 2023

Licensed under the MIT License, (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at:

    https://mit-license.org/

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

#import required packages
#importing packages to be used with error handling
try:
    import pandas as pd
except:
    print("Required package pandas not available")
    exit()
try:
    import numpy as np
except:
    print("Required package numpy not available")
    exit()
try:
    import argparse
except:
    print("Required package argparse not available")
    exit()
try:
    import subprocess
except:
    print("Required package subprocess not available")
    exit()
try:
    import os
except:
    print("Required package os not available")
    exit()

VIRGO2DBLoc = '/local/scratch/mfrance/VIRGO2/12_package/VIRGO2'

###Defining function which corrects mapping output based on percent of gene covered
def covCorr(samFile,coverageFile):

    coverage = pd.read_csv(coverageFile, sep='\t', header=0, usecols=[0, 5], names=['Gene','PercentCovered'])
    coverage = dict(zip(coverage['Gene'], coverage['PercentCovered']))

    #open the .sam file as a pandas dataframe
    sam = pd.read_csv(samFile, sep='\t', header=None, usecols=[0,1], names=['Read','Gene'],quotechar='@')
    sam = sam.groupby('Read')['Gene'].apply(list).to_dict()

    #a dictionary that keeps track of best gene for each read
    best_gene_per_read = {}

    for Read, Genes in sam.items():
        if len(Genes) == 1 and Genes[0] == '*':
            continue
        best_gene_per_read[Read] = max(Genes, key=lambda x : coverage[x])

    samOut = pd.DataFrame(data={'Read':list(best_gene_per_read.keys()), 'Gene':list(best_gene_per_read.values())})

    return samOut

def noCovCorr(samFile):
    samOut = pd.read_csv(samFile, delimiter='\t', header=None, usecols=[0,1], names=['Read','Gene'], low_memory=False)
    samOut = samOut[samOut['Gene']!='*']

    return samOut

def geneCounts(samOut):

    samOut=samOut.drop(['Read'],axis=1)
    samOut['Count']=1
    samOut = samOut.groupby("Gene").sum()
    samOut = samOut.sort_values(by='Count',ascending=False)
    #output to file
    return samOut

def garbageCollect(file2del):

    try:
        os.remove(file2del)
    except OSError as e:
        #Report failure to user
        print("Expected temporary files not available to delete\nExpected Files: {}".format(files2del))


####Argument Parsing
parser = argparse.ArgumentParser(description="VIRGO2 is a tool and associated database used to analyze vaginal shotgun metagenomes and metatranscriptomes")
subparsers = parser.add_subparsers(dest='command')

###Viable commands
#map -> map reads to VIRGO2 using bowtie2 and convert to gene counts
#compile -> compile mapping outputs from a directory of mapping results
#license -> print license

subparsers.required = True
#  subparser for mapping routine
parser_map = subparsers.add_parser('map')

# adding arguments for mapping routine
parser_map.add_argument('-r','--reads',help='Single-End reads file, can be gzipped',required=True)
parser_map.add_argument('-c','--cov',choices=['0','1'],help='Assign multi-mapped reads to gene with highest percent covered, 0=No,1=Yes, default:Yes',default='1',required=False)
parser_map.add_argument('-p','--threads',help='Number of threads used in mapping default:1',default=1,required=False)
parser_map.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the mapping output',required=True)
parser_map.add_argument('-b','--bypass',choices=['0','1'],help='Sam and coverager files already generated, proceed to coverage correction directly 0=No, 1=Yes, default No',default='0',required=False)
#subparser for compile routine
parser_compile = subparsers.add_parser('compile')
# adding arguments for compile routine
parser_compile.add_argument('-i','--input',help='Directory where bowtie2 mapping results are located')
parser_compile.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the compiled output')

#subparser for license routine
parser_license = subparsers.add_parser('license')

args = parser.parse_args()

#call functions for mapping
if args.command == 'map':

    if args.cov=='1':
        #####Mapping with coverage correction, DEFAULT
        if args.bypass=='0':
            try:
                subprocess.run(["bowtie2","-N","0","-L","20","-D","20","-R","3","-k","10","--sam-no-qname-trunc","--end-to-end","--seed","343",
                        "--no-unal","-p",args.threads,"-x",str(VIRGO2DBLoc),"-U",args.reads,"-S","{}.sam".format(args.outputPrefix)])
            except subprocess.CalledProcessError as e:
                print(e.output)

            #sorting sam file
            try:
                subprocess.run(['samtools','sort','{}.sam'.format(args.outputPrefix),'-o','{}.sorted.sam'.format(args.outputPrefix)])
            except subprocess.CalledProcessError as e:
                print(e.output)
        #generating gene level coverage values
            garbageCollect('{}.sam'.format(args.outputPrefix))

            try:
                subprocess.run(['samtools','coverage',"{}.sorted.sam".format(args.outputPrefix),'-o',"{}.cov".format(args.outputPrefix)])
            except subprocess.CalledProcessError as e:
                print(e.output)

        #####Viewing samfile to remove header
            try:
                #cut command prevents samtools from generating uneven number of columns for downstream read in pandas dataframe
                subprocess.run('samtools view -S {}.sorted.sam | cut -f 1,3 > {}.noHead.sam'.format(args.outputPrefix, args.outputPrefix), shell=True)
            
                #subprocess.run(['samtools','view',"{}.sorted.sam".format(args.outputPrefix),'-o',"{}.noHead.sam".format(args.outputPrefix)])
            except subprocess.CalledProcessError as e:
                print(e.output)
            garbageCollect('{}.sorted.sam'.format(args.outputPrefix))

        samOut = covCorr("{}.noHead.sam".format(args.outputPrefix),"{}.cov".format(args.outputPrefix))
        geneCounts = geneCounts(samOut)

        geneCounts.to_csv("{}.out".format(args.outputPrefix),sep="\t")

        #####Removing temporary files
        files2del=['{}.sam'.format(args.outputPrefix),'{}.sorted.sam'.format(args.outputPrefix),"{}.cov".format(args.outputPrefix),"{}.noHead.sam".format(args.outputPrefix)]
        garbageCollect('{}.noHead.sam'.format(args.outputPrefix))
        garbageCollect('{}.cov'.format(args.outputPrefix))
    elif args.cov=='0':
        #mapping without coverage correction, WARNING inflates gene counts
        try:
            subprocess.run(["bowtie2","-N","0","-L","20","-D","20","-R","3","--sam-no-qname-trunc","--end-to-end","--seed","343",
                        "--no-unal","-p",args.threads,"-x",str(VIRGO2DBLoc),"-U",args.reads,"-S","{}.sam".format(args.outputPrefix)])
        except subprocess.CalledProcessError as e:
            print(e.output)

        #sorting sam file
        try:
            subprocess.run(['samtools','sort','{}.sam'.format(args.outputPrefix),'-o','{}.sorted.sam'.format(args.outputPrefix)])
        except subprocess.CalledProcessError as e:
            print(e.output)
        garbageCollect('{}.sam'.format(args.outputPrefix))
        #####Viewing samfile to remove header
        try:
            subprocess.run(['samtools','view',"{}.sorted.sam".format(args.outputPrefix),'-o',"{}.noHead.sam".format(args.outputPrefix)])

        except subprocess.CalledProcessError as e:
            print(e.output)
        garbageCollect('{}.sorted.sam'.format(args.outputPrefix))
        samOut = noCovCorr("{}.noHead.sam".format(args.outputPrefix))
        geneCounts = geneCounts(samOut)
        garbageCollect('{}.noHead.sam'.format(args.outputPrefix))
        geneCounts.to_csv("{}.out".format(args.outputPrefix),sep="\t")

        files2del=['{}.sam'.format(args.outputPrefix),'{}.sorted.sam'.format(args.outputPrefix),"{}.noHead.sam".format(args.outputPrefix)]

#calls functions for compiling mapping results
if args.command == 'compile':

    output2concat = []
    for outFile in sorted(os.listdir(args.input)):
        print(outFile)
        try:
            if outFile.split('.')[-1]=='out':
                mappingDF = pd.read_csv(outFile,sep='\t',index_col=0)
                mappingDF.columns = [outFile.split('.')[0]]
                output2concat.append(mappingDF)
        except:
            continue
    try:
        compiledOutput = pd.concat(output2concat,axis=1)
        compiledOutput = compiledOutput.fillna(0)
        compiledOutput.to_csv("{}.summary.NR.txt".format(args.outputPrefix),sep="\t")
    except:
        print('No VIRGO2 mapping output files to compile in the designated directory')

if args.command == 'license':
    subprocess.run(['cat','License.txt'])
