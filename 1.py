#!/usr/bin/python3

import os , math , gzip , warnings , functools , random , fractions , numpy , sklearn.cluster , Bio.PDB , matplotlib.pyplot , mpl_toolkits.mplot3d
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

'''
#Identify number of residues in each secondary structure and their average phi and average psi angles
SecStruct = ''.join(SS)
#Loop
Loop = list()
for numb in SecStruct.replace('H' , '.').replace('S' , '.').split('.'):
	value = len(numb)
	if value == 0:
		pass
	else:
		Loop.append(value)
LoopNumb = len(Loop)
#PHI loop angels
PHI = [0]
angel = 0
for res , ang in zip(SS , Tor):			#Loop both lists
	if (res == 'H' or res == 'S'):		#If there is S or H from SS list ---> angel will become 0
		if angel != 0:			#And if the angel is a number other than 0
			PHI.append(angel)	#Append it
			angel = 0		#Make the angel 0
	else:					#If there is L from SS list
		angel = angel + ang[0]		#Make angel value from Tor[0] list
PHI.append(Tor[-1][0])
#PSI loop angels
PSI = list()
angel = 0
for res , ang in zip(SS , Tor):
	if (res == 'H' or res == 'S'):
		if angel != 0:
			PSI.append(angel)
			angel = 0
	else:
		angel = angel + ang[1]
PSI.append(0)
#Average loop angels
LOOP = list()
for L , Ph , Ps in zip(Loop , PHI , PSI):
	LOOP.append(((L) , (Ph/L) , (Ps/L)))
#Helix
helix = list()
for numb in SecStruct.replace('L' , '.').replace('S' , '.').split('.'):
	value = len(numb)
	if value == 0:
		pass
	else:
		helix.append(value)
helixNumb = len(helix)
#PHI helix angels
PHI = list()
angel = 0
for res , ang in zip(SS , Tor):
	if (res == 'L' or res == 'S'):
		if angel != 0:
			PHI.append(angel)
			angel = 0
	else:
		angel = angel + ang[0]
#PSI helix angels
PSI = list()
angel = 0
for res , ang in zip(SS , Tor):
	if (res == 'L' or res == 'S'):
		if angel != 0:
			PSI.append(angel)
			angel = 0
	else:
		angel = angel + ang[1]
#Average helix angels
HELIX = list()
for L , Ph , Ps in zip(helix , PHI , PSI):
	HELIX.append(((L) , (Ph/L) , (Ps/L)))
#Strand
strand = list()
for numb in SecStruct.replace('L' , '.').replace('H' , '.').split('.'):
	value = len(numb)
	if value == 0:
		pass
	else:
		strand.append(value)
strandNumb = len(strand)
#PHI strand angels
PHI = list()
angel = 0
for res , ang in zip(SS , Tor):
	if (res == 'L' or res == 'H'):
		if angel != 0:
			PHI.append(angel)
			angel = 0
	else:
		angel = angel + ang[0]
#PSI strand angels
PSI = list()
angel = 0
for res , ang in zip(SS , Tor):
	if (res == 'L' or res == 'H'):
		if angel != 0:
			PSI.append(angel)
			angel = 0
	else:
		angel = angel + ang[1]
#Average strand angels
STRAND = list()
for L , Ph , Ps in zip(strand , PHI , PSI):
	STRAND.append(((L) , (Ph/L) , (Ps/L)))





print('[+] GOOD\t' , TheFile)



'''
def DBSCAN(TheFile):
	''' Calculates the cluster centers of the Ramachandran plot, [output of the RamaPlot() function] '''
	''' Returns cluster centers '''
	df = pandas.read_csv('RamaPlot.csv' , sep=';')

	ML = sklearn.cluster.DBSCAN(eps = 0.8 , min_samples = 19)
	print(ML)


#--------------------------------------------------------------------------------------------------------------------------------------
DBSCAN('RamaPlot.csv')


'''
* DBSCAN center of clusters
* average angles plot
'''
