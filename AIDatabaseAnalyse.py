#!/usr/bin/python3

import os , math , gzip , Bio.PDB
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

def Database():
	''' A small script that cleans the PDB database, then isolates the secondary structure and the Phi/Psi torsion angles from each .pdb file '''
	''' Will generate the PDBDatabase directory with all the cleaned .pdb structures inside it, and the Data directory that contains the .csv files for all .pdb files '''
	From = int(smaller)
	To = int(bigger)
	#Collect Structures
	os.system('wget -rA .ent.gz ftp://ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/ -P DATABASE')
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
	#Separate Chains
	pdbfilelist = os.listdir('PDBDatabase')
	os.chdir('PDBDatabase')
	for thefile in pdbfilelist:
		#Open File
		TheFile = current + '/PDBDatabase/' + thefile
		TheName = TheFile.split('.')
		#Extract Each Chain and Save as Different Files
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
	#Remove Unwanted Structures
	pdbfilelist = os.listdir('PDBDatabase')
	ProteinCount = 1
	thedatafile = open('data.csv' , 'a')
	thedatafile.write(';1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;Distance_1;Distance_2;Distance_3;Distance_4;Distance_5;Distance_6;Distance_7;Distance_8;Distance_9;Distance_10;\n')
	current = os.getcwd()
	pdbfilelist = os.listdir('PDBDatabase')
	for thefile in pdbfilelist:
		TheFile = current + '/PDBDatabase/' + thefile
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only=False)
		#Delete Non-Protein Files
		if Type == []:
			print('[-] NOT PROTEIN\t' , thefile)
			os.remove(TheFile)
		else:
			#Get Secondary Structures
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
			else:
				#Renumber Residues
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
					final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line of the motif to have its atoms and residues sequencially labeled, as well as being in chain A
					AA2 = AA1
					PDB.write(final_line)													#Write to new file called motif.pdb
				PDB.close()
				os.remove(TheFile)
				os.rename(TheFile + 'X' , TheFile)
				#Get Torsion Angles
				count += 1
				Tor = list()
				for model in Bio.PDB.PDBParser().get_structure('X' , TheFile):
					for chain in model:
						polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
						for poly_index , poly in enumerate(polypeptides):
							phi_psi = poly.get_phi_psi_list()
							for res_index , residue in enumerate(poly):
								#Phi Angles
								if phi_psi[res_index][0] is None:
									phi = 0
								else:
									angle = phi_psi[res_index][0] * 180 / math.pi
									while angle > 180:
										angle = angle - 360
									while angle < -180:
										angle = angle + 360
									phi = angle
								#Psi Angles
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
				#Put Info Together
				ss = list()
				for val in SS:
					if val == 'L':
						ss.append('1')
					elif val == 'H':
						ss.append('2')
					elif val == 'S':
						ss.append('3')
				SecondaryStructures = ';' + ';'.join(ss)					#Secondary Structures L = 1, H = 2, S = 3 printed horisantally
				phiang = list()
				psiang = list()
				for val in Tor:
					phiang.append(val[0])
					psiang.append(val[1])
				PHIAngles = ';' + ';'.join(map(str, phiang))					#PHI angles printed horisantally
				PSIAngles = ';' + ';'.join(map(str, psiang))					#PSI angles printed horisantally
				Distances = ';' + ';'.join(map(str , distances))
				#Write To File
				line = str(ProteinCount) + SecondaryStructures + filling + Distances + '\n'	#The PHI and PSI angels are not being used because we cannot insert the angels as a feature during Machine Learning prediction, to use add this to the line variable: PHIAngles + filling + PSIAngles + filling
				thedatafile.write(line)
				ProteinCount += 1
				print('[+] GOOD\t' , thefile)
	thedatafile.close()
	os.system('rm -r PDBDatabase')





#--------------------------------------------------------------------------------------------------------------------------------------
Database
