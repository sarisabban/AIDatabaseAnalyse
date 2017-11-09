#!/usr/bin/python3

import os , gzip , warnings , Bio.PDB , matplotlib.pyplot , mpl_toolkits.mplot3d
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

def Database():
	''' This function downloads the full PDB database and cleans it up '''
	''' Out put will be a directory called PDBDatabase '''
	#Collect structures
	os.system('rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ ./DATABASE')
	current = os.getcwd()
	os.mkdir('PDBDatabase')
	filelist = os.listdir('DATABASE')
	for directories in filelist:
		files = os.listdir(current + '/DATABASE/' + directories)
		for afile in files:
			location = (current + '/DATABASE/' + directories + '/' + afile)
			print(location)
			os.rename(location , current + '/PDBDatabase/' + afile)
	os.system('rm -r ./DATABASE')
	#Clean Database
	pdbfilelist = os.listdir('PDBDatabase')
	io = Bio.PDB.PDBIO()
	os.chdir('PDBDatabase')
	for thefile in pdbfilelist:
		try:
			#Open file
			TheFile = current + '/PDBDatabase/' + thefile
			TheName = thefile.split('.')[0].split('pdb')[1].upper()
			#Extract file
			InFile = gzip.open(TheFile, 'rt')
			#Separate chains and save to different files
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure(TheName , InFile)
			count = 0
			for chain in structure.get_chains():
				io.set_structure(chain)
				io.save(structure.get_id() + '_' + chain.get_id() + '.pdb')
			print('[+] Extracted' + '\t' + thefile.upper())
			os.remove(TheFile)

		except:
			print('[-] Failed to Extracted' + '\t' + thefile.upper())
			os.remove(TheFile)
	os.chdir(current)
	#Remove unwanted structures
	current = os.getcwd()
	pdbfilelist = os.listdir('PDBDatabase')
	for thefile in pdbfilelist:
		TheFile = current + '/PDBDatabase/' + thefile
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure(TheFile.split('.')[0] , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only=True)
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

def RamaPlot(directory_to_search , plot_bool):
	''' Extract all the Phi and Psi angles of each residue in all structures in a directory and then plot them '''
	''' Generates the RamaPlot.pdf file '''
	current = os.getcwd()
	os.chdir(directory_to_search)
	Loop_phi = list()
	Loop_psi = list()
	Helix_phi = list()
	Helix_psi = list()
	Strand_phi = list()
	Strand_psi = list()
	#Identify the Phi and Psi angels for each secondary structure
	with warnings.catch_warnings(record=True) as w:							#Supress Bio.PDB user warnings
		for TheFile in os.listdir('.'):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure(TheFile.split('.')[0] , TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				for res in dssp:
					if res[2] == '-' or res[2] == 'T' or res[2] == 'S':		#Loop (DSSP code is - or T or S)
						Loop_phi.append(res[4])					#(PHI , PSI)
						Loop_psi.append(res[5])					#(PHI , PSI)
					elif res[2] == 'G' or res[2] == 'H' or res[2] == 'I':		#Helix (DSSP code is G or H or I)
						Helix_phi.append(res[4])				#(PHI , PSI)
						Helix_psi.append(res[5])				#(PHI , PSI)
					elif res[2] == 'B' or res[2] == 'E':				#Strand (DSSP code is B or E)
						Strand_phi.append(res[4])				#(PHI , PSI)
						Strand_psi.append(res[5])				#(PHI , PSI)
			except Exception as TheError:
				print(TheFile , '<---' , TheError)
			print('Ramachandran Plot - Got Phi and Psi Angels for\t' , TheFile.split('.')[0])
	#Replace 360 degrees with 0 degrees
	Loop_phi = [0.0 if x == 360.0 else x for x in Loop_phi]
	Loop_psi = [0.0 if x == 360.0 else x for x in Loop_psi]
	Helix_phi = [0.0 if x == 360.0 else x for x in Helix_phi]
	Helix_psi = [0.0 if x == 360.0 else x for x in Helix_psi]
	Strand_phi = [0.0 if x == 360.0 else x for x in Strand_phi]
	Strand_psi = [0.0 if x == 360.0 else x for x in Strand_psi]
	#Plot full graph
	if plot_bool == 1:
		matplotlib.rcParams['axes.facecolor'] = '0.8'
		matplotlib.pyplot.scatter(Loop_phi, Loop_psi , label='Loop' , color = '#038125' , s = 0.5)
		matplotlib.pyplot.scatter(Helix_phi, Helix_psi , label='Helix' , color = '#c82100' , s = 0.5)
		matplotlib.pyplot.scatter(Strand_phi, Strand_psi , label='Strand' , color = '#ffe033' , s = 0.5)
		matplotlib.pyplot.legend(bbox_to_anchor = (0. , 1.02 , 1. , .102) , loc = 3 , ncol = 3 , mode = 'expand' , borderaxespad = 0. , facecolor = '0.9' , markerscale = 10)
		matplotlib.pyplot.title('Ramachandran Plot' , y = 1.08)
		matplotlib.pyplot.xlabel('Phi Angels')
		matplotlib.pyplot.ylabel('Psi Angels')
		matplotlib.pyplot.ylim(-180 , 180)
		matplotlib.pyplot.xlim(-180 , 180)
		os.chdir(current)
		matplotlib.pyplot.savefig('RamaPlot.pdf')
	else:
		pass

def Numbers(directory_to_search , plot_bool , show_plot):
	''' Count the number of secondary structures in each protein '''
	''' Returns a tuple of lists with loop numbers [0], helix number [1], and strand number [2] '''
	current = os.getcwd()
	os.chdir(directory_to_search)
	LoopLen = list()
	HelixLen = list()
	StrandLen = list()
	with warnings.catch_warnings(record=True) as w:							#Supress Bio.PDB user warnings
		for TheFile in os.listdir('.'):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure(TheFile.split('.')[0] , TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				Loop = list()
				Helix = list()
				Strand = list()
				for res in dssp:
					if res[2] == '-' or res[2] == 'T' or res[2] == 'S':		#Loop (DSSP code is - or T or S)
						Loop.append('L')
					elif res[2] == 'G' or res[2] == 'H' or res[2] == 'I':		#Helix (DSSP code is G or H or I)
						Helix.append('H')
					elif res[2] == 'B' or res[2] == 'E':				#Strand (DSSP code is B or E)
						Strand.append('S')
				LoopLen.append(len(Loop))
				HelixLen.append(len(Helix))
				StrandLen.append(len(Strand))
				Loop = list()
				Helix = list()
				Strand = list()
			except Exception as TheError:
				print(TheFile , '<---' , TheError)
			print('Counted Secondary Structures for\t' , TheFile.split('.')[0])
	os.chdir(current)
	#Plot graph
	if plot_bool == 1:
		matplotlib.rcParams['axes.facecolor'] = '0.5'
		AX = mpl_toolkits.mplot3d.Axes3D(matplotlib.pyplot.figure())
		AX.scatter(LoopLen , HelixLen , StrandLen , c = 'red' , s = 0.5)
		AX.set_title('Number of Secondary Structures Plot' , y = 1.02)
		AX.set_xlabel('Number of Loops')
		AX.set_ylabel('Number of Helices')
		AX.set_zlabel('Number of Strands')
		if show_plot == 1:
			matplotlib.pyplot.savefig('SSnumPlot.pdf')
			matplotlib.pyplot.show()
		else:
			matplotlib.pyplot.savefig('SSnumPlot.pdf')
	else:
		pass
	return(LoopLen , HelixLen , StrandLen)
#--------------------------------------------------------------------------------------------------------------------------------------
#Database()
#RamaPlot('PDBDatabase' , 1)
