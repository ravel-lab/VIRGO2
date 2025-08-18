#!/usr/bin/python
import pandas as pd 
import sys
import numpy as np
from datetime import date
import os


#this script takes as input the summary.NR output file from VIRGO
#adds to it the Taxa, KEGG, and Length information by merging with
#VIRGO annnotation tables

#collecting the date for use in output filename
today = date.today().strftime("%m%d%y")

#file prefix as second sys arguement
prefix = str(sys.argv[2])

data = pd.read_csv(sys.argv[1],sep="\t")
#data.columns = ['PC']
data['Gene'] = data['Gene'].astype('str')

scriptLoc = os.path.realpath(__file__)
scriptLoc = scriptLoc.split('VIRGO2_add_annotations.py')[0]


#gene length data, every gene has a length
len_data = pd.read_csv("{}../AnnotationTables/0.VIRGO2.geneLength.txt".format(scriptLoc),sep="\t",header=None)
len_data.columns = ['Gene','Length']
len_data['Gene'] = len_data['Gene'].astype('str')

#taxa annotation, not every gene has designation
taxa_data = pd.read_csv("{}../AnnotationTables/1.VIRGO2.taxon.txt".format(scriptLoc),sep="\t")
taxa_data = taxa_data.drop(['Cluster','Cat'],axis=1)
taxa_data['Gene'] = taxa_data['Gene'].astype('str')

#kegg annotation, not every gene has a designation
kegg_data = pd.read_csv("{}../AnnotationTables/7.VIRGO2.kegg.txt".format(scriptLoc),sep="\t")
kegg_data = kegg_data[['Gene','KEGG']]
kegg_data['Gene'] = kegg_data['Gene'].astype('str')

#geneProduct annotation, not every gene has annotation
gp_data = pd.read_csv("{}../AnnotationTables/6.VIRGO2.geneProduct.txt".format(scriptLoc),sep="\t")
gp_data['Gene'] = gp_data['Gene'].astype('str')


data = pd.merge(left=len_data,right=data,left_on="Gene",right_on='Gene',how='right')
data = pd.merge(left=taxa_data,right=data,left_on="Gene",right_on='Gene',how='right')
data = pd.merge(left=kegg_data,right=data,on="Gene",how='right')
data = pd.merge(left=gp_data,right=data,on="Gene",how='right')

data.to_csv("%s_NR_anno_%s.csv" %(prefix,today),sep=",",index=None)