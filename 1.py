#!/usr/bin/python3

import os , math , gzip , warnings , random , Bio.PDB , matplotlib.pyplot
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

def Dataset():
	''' Loops through all structures of the database and extracts information to generate a dataset for Machine Learning '''
	''' Generates the Data.set file '''
	current = os.getcwd()
	DATA = open('data.csv' , 'a')
	DATA.write(';Secondary Structure;Label;Number of Residues;Average Phi;Average Psi\n')
	os.chdir('PDBDatabase')
	count = 0
	with warnings.catch_warnings(record=True) as w:					#Supress Bio.PDB user warnings
		for TheFile in os.listdir('.'):
			try:
				#Identify secondary structures
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure('X' , TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				SS = list()
				for res in dssp:
					ss = res[2]
					if ss == '-' or ss == 'T' or ss == 'S':		#Loop (DSSP code is - or T or S)
						SS.append('L')
					elif ss == 'G' or ss == 'H' or ss == 'I':	#Helix (DSSP code is G or H or I)
						SS.append('H')
					elif ss == 'B' or ss == 'E':			#Sheet (DSSP code is B or E)
						SS.append('S')
			except Exception as TheError:
				print(TheFile , '<---------------' , TheError)
			#Get torsion angles
			Tor = list()
			for model in Bio.PDB.PDBParser().get_structure('X' , TheFile):
				for chain in model:
					polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
					for poly_index , poly in enumerate(polypeptides):
						phi_psi = poly.get_phi_psi_list()
						for res_index , residue in enumerate(poly):
							#Phi angles
							if phi_psi[res_index][0] is None:
								phi = 0
							else:
								angle = phi_psi[res_index][0] * 180 / math.pi
								while angle > 180:
									angle = angle - 360
								while angle < -180:
									angle = angle + 360
								phi = angle
							#Psi angles
							if phi_psi[res_index][1] is None:
								psi = 0
							else:
								angle = phi_psi[res_index][1] * 180 / math.pi
								while angle > 180:
									angle = angle - 360
								while angle < -180:
									angle = angle + 360
								psi = angle
							Tor.append((phi , psi))
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
			#Generate Data.set file
			#Add loop info
			Lines = list()
			for item in LOOP:
				L = ((';' + 'Loop' + ';' + '0' + ';' + str(item[0]) + ';' + str(item[1]) + ';' + str(item[2]) + '\n'))
				Lines.append(L)
			#Add helix info
			for item in HELIX:
				H = ((';' + 'Helix' + ';' + '1' + ';' + str(item[0]) + ';' + str(item[1]) + ';' + str(item[2]) + '\n'))
				Lines.append(H)
			#Add strand info
			for item in STRAND:
				S = ((';' + 'Strand' + ';' + '2' + ';' + str(item[0]) + ';' + str(item[1]) + ';' + str(item[2]) + '\n'))
				Lines.append(S)
			#Randomise entries and write to file
			random.shuffle(Lines)
			for line in Lines:
				count += 1
				DATA.write((str(count) + line))
			print('[+] GOOD\t' , TheFile)
	DATA.close()
	os.chdir(current)

























































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












#--------------------------------------------------------------------------------------------------------------------------------------
RamaPlot('PDBDatabase' , 1)



'''
2. average angles plot
3. angles of each secondary structure: purple S   red H    green L
4. number of secondary structures in each protein
5. ratio of helix to strand to loop to each other and to size of protein
6. leangth of helix and strand to size of protein
'''
