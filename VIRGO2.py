#######VIRGO2
#######Authors: Michael France, Issac Chaudry
#######Contact: mfrance@som.umaryland.edu

####DEPENDENCIES required to be in PATH
###bowtie2
###samtools

###V1 081825

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

scriptLoc = os.path.realpath(__file__)
scriptLoc = scriptLoc.split('VIRGO2.py')[0]

VIRGO2DB = '{}/Index/VIRGO2'.format(scriptLoc)

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
#install -> unzip annotation and fasta files and build bowtie2 index
#map -> map reads to VIRGO2 using bowtie2 and convert to gene counts
#compile -> compile mapping outputs from a directory of mapping results
#taxonomy -> produce taxonomic composition estimates from compiled VIRGO2 output
#license -> print license

subparsers.required = True

#subparser for install command
parser_install = subparsers.add_parser('install')
parser_install.add_argument('-p','--threads',help='Number of threads used in building index default:1',default='1',required=False)


#subparser for mapping routine
parser_map = subparsers.add_parser('map')

# adding arguments for mapping routine
parser_map.add_argument('-r','--reads',help='Single-End reads file, can be gzipped',required=True)
parser_map.add_argument('-c','--cov',choices=['0','1'],help='Assign multi-mapped reads to gene with highest percent covered, 0=No,1=Yes, default:Yes',default='1',required=False)
parser_map.add_argument('-p','--threads',help='Number of threads used in mapping default:1',default='1',required=False)
parser_map.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the mapping output',required=True)
parser_map.add_argument('-b','--bypass',choices=['0','1'],help='Sam and coverage files already generated, proceed to coverage correction directly 0=No, 1=Yes, default No',default='0',required=False)
#subparser for compile routine
parser_compile = subparsers.add_parser('compile')
# adding arguments for compile routine
parser_compile.add_argument('-i','--input',help='Full path to directory where bowtie2 mapping results are located')
parser_compile.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the compiled output')

parser_taxonomy = subparsers.add_parser('taxonomy')
parser_taxonomy.add_argument('-i','--input',help='Full path to compiled results file')
parser_taxonomy.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the compiled output')
parser_taxonomy.add_argument('-b','--bacteria',choices=['0','1'],help='Report composition including only the bacteria, default=1',default='1',required=False)
parser_taxonomy.add_argument('-f','--filter',choices=['0','1'],help='Mask contribution from taxa where number of genes detected is below threshold, default=1',default='1',required=False)
parser_taxonomy.add_argument('-r','--readCounts',choices=['0','1'],help='Report values as read counts instead of relative abundances, default=0',default='0',required=False)
parser_taxonomy.add_argument('-m','--multigenera',choices=['0','1'],help="Report relative abundance of multigenera genes, off by default (0)",default='0',required=False)

#subparser for license routine
parser_license = subparsers.add_parser('license')

args = parser.parse_args()

if args.command == 'install':

    for outFile in sorted(os.listdir("{}FastaFiles/".format(scriptLoc))):
        print(outFile)
        try:
            if outFile.split('.')[-1]=='gz':
                subprocess.run(["gunzip","{path}FastaFiles/{file}".format(path=scriptLoc,file=outFile)])
            else:
                pass    
        except subprocess.CalledProcessError as e:
            print(e.output)

    for outFile in sorted(os.listdir("{}AnnotationTables/".format(scriptLoc))):
        print(outFile)
        try:
            if outFile.split('.')[-1]=='gz':
                subprocess.run(["gunzip","{path}AnnotationTables/{file}".format(path=scriptLoc,file=outFile)])
            else:
                pass    
        except subprocess.CalledProcessError as e:
            print(e.output)

    try:
        subprocess.run(["bowtie2-build","-f","--seed","343","--threads",args.threads,"{path}FastaFiles/VIRGO2.fa".format(path=scriptLoc),"{}Index/VIRGO2".format(scriptLoc)])
    except subprocess.CalledProcessError as e:
        print(e.output)

    try:
        subprocess.run(["rm","{}Index/.placeholder".format(scriptLoc)])
    except subprocess.CalledProcessError as e:
        print(e.output)


#call functions for mapping
if args.command == 'map':

    if args.cov=='1':
        #####Mapping with coverage correction, DEFAULT
        if args.bypass=='0':
            try:
                subprocess.run(["bowtie2","-N","0","-L","20","-D","20","-R","3","-k","10","--sam-no-qname-trunc","--end-to-end","--seed","343",
                        "--no-unal","-p",args.threads,"-x",str(VIRGO2DB),"-U",args.reads,"-S","{}.sam".format(args.outputPrefix)])
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
                        "--no-unal","-p",args.threads,"-x",str(VIRGO2DB),"-U",args.reads,"-S","{}.sam".format(args.outputPrefix)])
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
                mappingDF = pd.read_csv("{path}/{file}".format(path=str(args.input),file=str(outFile)),sep='\t',index_col=0)
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

if args.command == 'taxonomy':

    print(args.input)
    
    #####reading in compiled output
    try:
        VIRGO2table = pd.read_csv(args.input,sep='\t')
        VIRGO2table['Gene']=VIRGO2table['Gene'].astype(str)
    except:
        print('No compiled VIRGO2 output provided')

    #####reading in annotation tables
    try:
        lenKey = pd.read_csv("{}/AnnotationTables/0.VIRGO2.geneLength.txt".format(scriptLoc),sep="\t")
        lenKey['Gene'] = lenKey['Gene'].astype('str')
    except:
        print("Could not read the gene length annotation table (0.VIRGO2.geneLength.txt)")
    #taxa annotation, not every gene has designation
    try:
        taxKey = pd.read_csv("{}/AnnotationTables/1.VIRGO2.taxon.txt".format(scriptLoc),sep="\t")
        taxKey['Gene'] = taxKey['Gene'].astype('str')
        taxKey = taxKey.drop(['Cluster','Cat'],axis=1)
    except:
        print("Could not read the taxon annotation table (1.VIRGO2.taxon.txt)")

    VIRGO2table = pd.merge(left=lenKey,right=VIRGO2table,left_on="Gene",right_on='Gene',how='right')
    VIRGO2table = pd.merge(left=taxKey,right=VIRGO2table,left_on="Gene",right_on='Gene',how='right')

    try:
        thresh = pd.read_csv("{}/AnnotationTables/2.VIRGO2.taxonThresholds.txt".format(scriptLoc),sep="\t",index_col=0)
        bacteriaList=list(thresh[thresh['Domain']=='Bacteria'].index)
        thresh = thresh.drop(['Domain','MedianGenes'],axis=1)
        thresh = thresh[thresh['Thresh']>0]
    except:
        print("Could not read taxon threshold file (2.VIRGO2.taxonThresholds.txt)")
    
    #removes non-bacterial taxa by default
    if args.bacteria=='1':
        VIRGO2table=VIRGO2table[VIRGO2table['Taxa'].isin(bacteriaList)]
    elif args.bacteria=='0':
        pass


    #counting number of genes detected per species per metagenome
    VIRGO2bool = pd.concat([VIRGO2table[VIRGO2table.columns[1:2]],VIRGO2table[VIRGO2table.columns[3:]].astype(bool)],axis=1).groupby(by='Taxa').sum()
    VIRGO2bool = pd.merge(left=thresh,right=VIRGO2bool,left_index=True,right_index=True,how="inner")
    VIRGO2bool = VIRGO2bool[VIRGO2bool.columns[1:]].ge(VIRGO2bool['Thresh'],axis=0)

    VIRGO2table = VIRGO2table.drop(['Gene'],axis=1)
    

    if args.multigenera=='0':
        VIRGO2table = VIRGO2table[VIRGO2table['Taxa']!='MultiGenera']
        VIRGO2bool=VIRGO2bool.drop(['MultiGenera'],axis=0)
    elif args.multigenera=='1':
        pass

    VIRGO2table = VIRGO2table[VIRGO2table['Taxa'].notna()]

    VIRGO2table = VIRGO2table.groupby("Taxa").sum()
    VIRGO2table = VIRGO2table.reset_index()
    VIRGO2table = VIRGO2table[VIRGO2table['Taxa'].isin(VIRGO2bool.index)]
    VIRGO2table = VIRGO2table.drop(['Length'],axis=1)

    #filtering data by detected taxa if requested, happens by default
    if args.filter=='1':
        VIRGO2table = VIRGO2table[VIRGO2table['Taxa'].isin(VIRGO2bool.index)]
        VIRGO2table = VIRGO2table.set_index("Taxa")
        TaxaTableOut = VIRGO2table*VIRGO2bool
    elif args.filter=='0':
        TaxaTableOut = VIRGO2table.set_index("Taxa")
    #calculating relative abundance, if requested, happens by default

    if args.readCounts=='0':
        TaxaTableOut=TaxaTableOut.div(TaxaTableOut.sum(axis=0),axis=1)
        reportType='relAbund'
    elif args.readCounts=='1':
        TaxaTableOut=TaxaTableOut
        reportType='readCount'

    #sorting the output taxa by average relative abundance across all samples
    TaxaTableOut['taxa_ave_abundance'] = TaxaTableOut.sum(axis=1)
    TaxaTableOut = TaxaTableOut.sort_values(by=["taxa_ave_abundance"],ascending=False)
    TaxaTableOut = TaxaTableOut.drop(['taxa_ave_abundance'],axis=1)
    
    TaxaTableOut = TaxaTableOut.T
    TaxaTableOut.to_csv("{}.{}.csv".format(args.outputPrefix,reportType),sep=",")

if args.command == 'license':
    subprocess.run(['cat','{}/LICENSE'.format(scriptLoc)])
