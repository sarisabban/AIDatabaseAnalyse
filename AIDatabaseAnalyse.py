#!/usr/bin/python3

import os , math , gzip , warnings , Bio.PDB
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
	os.chdir('PDBDatabase')
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

			print(SS)
			print(Tor)














		'''
			#Separate secondary structures
			SecStruct = ''.join(SS)
			Loop = list()
			for numb in SecStruct.replace('H' , '.').replace('S' , '.').split('.'):
				value = len(numb)
				if value == 0:
					pass
				else:
					Loop.append(value)
			LoopNumb = len(Loop)
			Helix = list()
			for numb in SecStruct.replace('L' , '.').replace('S' , '.').split('.'):
				value = len(numb)
				if value == 0:
					pass
				else:
					Helix.append(value)
			HelixNumb = len(Helix)
			Strand = list()
			for numb in SecStruct.replace('L' , '.').replace('H' , '.').split('.'):
				value = len(numb)
				if value == 0:
					pass
				else:
					Strand.append(value)
			StrandNumb = len(Strand)

			print(LoopNumb , Loop , HelixNumb , Helix, StrandNumb , Strand)
		'''











































def ML():
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
Dataset()
