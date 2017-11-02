#!/usr/bin/python3

import os , math , gzip , warnings , random , Bio.PDB
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

def Database():
	''' This function downloads the full PDB database and cleans it up '''
	''' Out put will be a directory called PDBDatabase '''
	#Collect structures
#	os.system('wget -rA .ent.gz ftp://ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/ -P DATABASE')
	current = os.getcwd()
	os.mkdir('PDBDatabase')
	filelist = os.listdir('DATABASE/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb')
	for directories in filelist:
		files = os.listdir(current + '/DATABASE/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/' + directories)
		for afile in files:
			location = (current + '/DATABASE/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/' + directories + '/' + afile)
			print(location)
			os.rename(location , current + '/PDBDatabase/' + afile)
	os.system('rm -r ./DATABASE')
	#Separate chains
	pdbfilelist = os.listdir('PDBDatabase')
	os.chdir('PDBDatabase')
	for thefile in pdbfilelist:
		#Open file
		TheFile = current + '/PDBDatabase/' + thefile
		TheName = TheFile.split('.')
		#Extract each chain and save as a different file
		InFile = gzip.open(TheFile, 'rb')
		for line in InFile:
			line = line.decode()
			if line.startswith('ATOM'):
				chain = line[21]
				Name = TheName[0].split('pdb')
				output = open(Name[1].upper() + '_' + chain + '.pdb' , 'a')
				output.write(line)
				output.close()
			else:
				pass
		print('[+] Extracted' + '\t' + Name[1].upper() , '\t' , chain)
		os.remove(TheFile)
	os.chdir(current)
	#Remove unwanted structures
	current = os.getcwd()
	pdbfilelist = os.listdir('PDBDatabase')
	for thefile in pdbfilelist:
		TheFile = current + '/PDBDatabase/' + thefile
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only=False)
		#Delete non-protein files
		if Type == []:
			print('[-] NOT PROTEIN\t' , thefile)
			os.remove(TheFile)
		else:
			#Renumber residues
			pdb = open(TheFile , 'r')
			PDB = open(TheFile + 'X' , 'w')
			count = 0
			num = 0
			AA2 = None
			for line in pdb:
				count += 1														#Sequencially number atoms
				AA1 = line[23:27]													#Sequencially number residues
				if not AA1 == AA2:
					num += 1			
				final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line to have its atoms and residues sequencially labeled, as well as being in chain A
				AA2 = AA1
				PDB.write(final_line)													#Write to new file called motif.pdb
			PDB.close()
			print('[+] GOOD\t' , thefile)
			os.remove(TheFile)
			os.rename(TheFile + 'X' , TheFile)

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

def RamaPLot():
	current = os.getcwd()
	DATA = open('Rama' , 'a')
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
			
							DATA.write((str(phi) + '\t' + str(psi) + '\n'))
			print('[+] GOOD\t' , TheFile)					
	DATA.close()
	os.chdir(current)






def GnuPlot():
	pass




def ML():
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
#Dataset()
#RamaPLot()

data = open('data.csv' , 'r')
files = open('data' , 'a')
for line in data:
	line = line.split(';')
	files.write((line[4] + '\t' + line[5]))
data.close()
files.close()

'''
set xrange [-180:180]
set yrange [-180:180]
plot 'Rama' with dots lc rgb 'red' , 'data' with dots lc rgb 'blue'
'''
