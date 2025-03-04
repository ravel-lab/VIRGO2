#!/usr/bin/python
import pandas as pd 
import sys
import numpy as np

#importing the datafile
data = pd.read_csv(sys.argv[1],sep=",")

#dividing each gene's read counts by its legnth, to correct
glcorr_data = data[data.columns[4:]].div(data['Length'],axis=0)

data = pd.concat([data[data.columns[0:4]],glcorr_data],axis=1)
data = data[data['Cat']!='MultiGenera']
data = data[data['Cat'].notna()]
data = data.drop(['Cat','Length','Gene','GeneProduct'],axis=1)

data_taxa = data.groupby('Taxa').sum()

data_taxa_rel = data_taxa.div(data_taxa.sum(axis=0),axis=1)

data_taxa_rel['taxa_abundance'] = data_taxa_rel.sum(axis=1)
data_taxa_rel = data_taxa_rel.sort_values(by=["taxa_abundance"],ascending=False)
data_taxa_rel = data_taxa_rel.drop(['taxa_abundance'],axis=1)
data_taxa_rel.to_csv(sys.argv[2],sep=",")
