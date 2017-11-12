#!/usr/bin/python3

import pandas , numpy , sklearn.cluster , matplotlib.pyplot
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

def DBSCAN(TheFile):
	''' Calculates the cluster centers of the Ramachandran plot, [output of the RamaPlot() function] '''
	''' Returns cluster centers '''
	df = pandas.read_csv(TheFile , sep=';')
	L = df[['Loop PHI' , 'Loop PSI']]
	H = df[['Helix PHI' , 'Helix PSI']]
	S = df[['Strand PHI' , 'Strand PSI']]
	ML = sklearn.cluster.DBSCAN(eps = 2 , min_samples = 10).fit(L)
	print(ML.labels_)


#--------------------------------------------------------------------------------------------------------------------------------------
#DBSCAN('RamaPlot.csv')
