#!/usr/bin/python3

import os , gzip , warnings , itertools , numpy , Bio.PDB , matplotlib.pyplot , mpl_toolkits.mplot3d , pandas
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
	''' Generates a Ramachandran plot called RamaPlot.pdf and saves the data for futher analysis Loop_phi [0] , Loop_psi [1] , Helix_phi [2] , Helix_psi [3] , Strand_phi [4] , and Strand_psi [5] in the RamaPlot.csv file'''
	data = open('RamaPlot.csv' , 'a')
	data.write('Loop PHI;Loop PSI;Helix PHI;Helix PSI;Strand PHI;Strand PSI\n')
	current = os.getcwd()
	os.chdir(directory_to_search)
	#Identify the Phi and Psi angels for each secondary structure
	with warnings.catch_warnings(record=True) as w:							#Supress Bio.PDB user warnings
		for TheFile in os.listdir('.'):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure(TheFile.split('.')[0] , TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				Loop_phi = list()
				Loop_psi = list()
				Helix_phi = list()
				Helix_psi = list()
				Strand_phi = list()
				Strand_psi = list()
				for res in dssp:
					if res[2] == '-' or res[2] == 'T' or res[2] == 'S':		#Loop (DSSP code is - or T or S)
						Loop_phi.append(res[4])					#(PHI)
						Loop_psi.append(res[5])					#(PSI)
					elif res[2] == 'G' or res[2] == 'H' or res[2] == 'I':		#Helix (DSSP code is G or H or I)
						Helix_phi.append(res[4])				#(PHI)
						Helix_psi.append(res[5])				#(PSI)
					elif res[2] == 'B' or res[2] == 'E':				#Strand (DSSP code is B or E)
						Strand_phi.append(res[4])				#(PHI)
						Strand_psi.append(res[5])				#(PSI)
				#Replace 360 degrees with 0 degrees
				Loop_phi = [0.0 if x == 360.0 else x for x in Loop_phi]
				Loop_psi = [0.0 if x == 360.0 else x for x in Loop_psi]
				Helix_phi = [0.0 if x == 360.0 else x for x in Helix_phi]
				Helix_psi = [0.0 if x == 360.0 else x for x in Helix_psi]
				Strand_phi = [0.0 if x == 360.0 else x for x in Strand_phi]
				Strand_psi = [0.0 if x == 360.0 else x for x in Strand_psi]
				for Lph , Lps , Hph , Hps , Sph , Sps in itertools.zip_longest(Loop_phi , Loop_psi , Helix_phi , Helix_psi , Strand_phi , Strand_psi , fillvalue = '0'):
					line = str(Lph) + ';' + str(Lps) + ';' + str(Hph) + ';' + str(Hps) + ';' + str(Sph) + ';' + str(Sps) + '\n'
					data.write(line)
			except Exception as TheError:
				print(TheFile , '<---' , TheError)
			print('Ramachandran plot - got phi and psi angels for\t' , TheFile.split('.')[0])
	data.close()
	os.chdir(current)
	#Plot full graph
	if plot_bool == 1:
		df = pandas.read_csv('RamaPlot.csv' , sep=';')
		matplotlib.rcParams['axes.facecolor'] = '0.8'
		L = df.plot(x = 'Loop PHI' , y = 'Loop PSI' , kind = 'scatter' , label = 'Loop' , color = '#038125' , s = 0.5)
		H = df.plot(x = 'Helix PHI' , y = 'Helix PSI' , kind = 'scatter'  , label = 'Helix' , color = '#c82100' , s = 0.5 , ax = L)
		S = df.plot(x = 'Strand PHI' , y = 'Strand PSI' , kind = 'scatter'  , label = 'Strand' , color = '#ffe033' , s = 0.5 , ax = H)
		matplotlib.pyplot.legend(bbox_to_anchor = (0. , 1.02 , 1. , .102) , loc = 3 , ncol = 3 , mode = 'expand' , borderaxespad = 0. , facecolor = '0.9' , markerscale = 10)
		matplotlib.pyplot.title('Ramachandran Plot' , y = 1.08)
		matplotlib.pyplot.xlabel('Phi Angels')
		matplotlib.pyplot.ylabel('Psi Angels')
		matplotlib.pyplot.ylim(-180 , 180)
		matplotlib.pyplot.xlim(-180 , 180)
		matplotlib.pyplot.savefig('RamaPlot.pdf')
	else:
		pass

def Numbers(directory_to_search , plot_bool , show_plot):
	''' Count the number of secondary structures in each protein '''
	''' Generates a plot with the number of each secondary structure in each protein called SSnumPlot and returns a tuple of lists with loop numbers [0], helix number [1], strand number [2], protein's size [3] '''
	current = os.getcwd()
	os.chdir(directory_to_search)
	LoopLen = list()
	HelixLen = list()
	StrandLen = list()
	Sizelen = list()
	with warnings.catch_warnings(record=True) as w:							#Supress Bio.PDB user warnings
		for TheFile in os.listdir('.'):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure(TheFile.split('.')[0] , TheFile)
				model = structure[0]
				ppb = Bio.PDB.Polypeptide.PPBuilder()
				Type = ppb.build_peptides(structure , aa_only=True)
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				for aa in dssp:
					length = aa[0]
				Sizelen.append(length)
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
			print('Counted secondary structures for\t' , TheFile.split('.')[0])
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
	return(LoopLen , HelixLen , StrandLen , Sizelen)

def Probability(TheList , plot_bool , show_plot):
	''' Calculates the probability of secondary structures given their total numbers and protein sizes [output from the Numbers() function] '''
	''' Generates a plot with the probability of each secondary structure called SSnumPlot.pdf '''
	loopfrac = list()
	helixfrac = list()
	strandfrac = list()	
	for loop , helix , strand , size in zip(TheList[0] , TheList[1] , TheList[2] , TheList[3]):
		loopfrac.append(float(loop / size))
		helixfrac.append(float(helix / size))
		strandfrac.append(float(strand / size))
	#Plot graph
	if plot_bool == 1:
		matplotlib.rcParams['axes.facecolor'] = '0.5'
		AX = mpl_toolkits.mplot3d.Axes3D(matplotlib.pyplot.figure())
		AX.scatter(loopfrac , helixfrac , strandfrac , c = 'red' , s = 0.5)
		AX.set_title('Probability of Secondary Structures Plot' , y = 1.02)
		AX.set_xlabel('Probability of Loops')
		AX.set_ylabel('Probability of Helices')
		AX.set_zlabel('Probability of Strands')
		if show_plot == 1:
			matplotlib.pyplot.savefig('SSproPlot.pdf')
			matplotlib.pyplot.show()
		else:
			matplotlib.pyplot.savefig('SSproPlot.pdf')
	else:
		pass

def Length(directory_to_search , plot_bool , show_plot):
	''' Caulculates the average length of secondary structures within each protein '''
	''' Generates a plot with the average length of each secondary structure in each protein called SSlenPlot.pdb '''
	current = os.getcwd()
	os.chdir(directory_to_search)
	Loopav = list()
	Helixav = list()
	Strandav = list()
	Sizenum = list()	
	with warnings.catch_warnings(record=True) as w:							#Supress Bio.PDB user warnings
		for TheFile in os.listdir('.'):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure(TheFile.split('.')[0] , TheFile)
				model = structure[0]
				ppb = Bio.PDB.Polypeptide.PPBuilder()
				Type = ppb.build_peptides(structure , aa_only=True)
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				SS = list()
				Size = None
				#Identify length of structure
				for aa in dssp:
					Size = aa[0]
				#Identify secondary structures
				for res in dssp:
					if res[2] == '-' or res[2] == 'T' or res[2] == 'S':		#Loop (DSSP code is - or T or S)
						SS.append('L')
					elif res[2] == 'G' or res[2] == 'H' or res[2] == 'I':		#Helix (DSSP code is G or H or I)
						SS.append('H')
					elif res[2] == 'B' or res[2] == 'E':				#Strand (DSSP code is B or E)
						SS.append('S')
				FinalSS = ''.join(SS)
				#Loop
				Loop = list()
				for numb in FinalSS.replace('H' , '.').replace('S' , '.').split('.'):
					value = len(numb)
					if value == 0:
						pass
					else:
						Loop.append(value)
				if Loop == []: Loop.append(0)
				#Helix
				Helix = list()
				for numb in FinalSS.replace('L' , '.').replace('S' , '.').split('.'):
					value = len(numb)
					if value == 0:
						pass
					else:
						Helix.append(value)
				if Helix == []: Helix.append(0)
				#Strand
				Strand = list()
				for numb in FinalSS.replace('H' , '.').replace('L' , '.').split('.'):
					value = len(numb)
					if value == 0:
						pass
					else:
						Strand.append(value)
				if Strand == []: Strand.append(0)
				Loopav.append(numpy.mean(Loop))
				Helixav.append(numpy.mean(Helix))
				Strandav.append(numpy.mean(Strand))
				Sizenum.append(Size)
			except Exception as TheError:
				print(TheFile , '<---' , TheError)
			print('Counted average length of all secondary structures for\t' , TheFile.split('.')[0])
	os.chdir(current)
	#Plot graph
	if plot_bool == 1:
		matplotlib.rcParams['axes.facecolor'] = '0.5'
		AX = mpl_toolkits.mplot3d.Axes3D(matplotlib.pyplot.figure())
		AX.scatter(Sizenum , Helixav , Strandav , c = 'red' , s = 0.5)
		AX.set_title('Number of Helices and Strands Plot' , y = 1.02)
		AX.set_xlabel('Size of Protein')
		AX.set_ylabel('Average Length of Helices')
		AX.set_zlabel('Average Length of Strands')
		if show_plot == 1:
			matplotlib.pyplot.savefig('SSlenPlot.pdf')
			matplotlib.pyplot.show()
		else:
			matplotlib.pyplot.savefig('SSlenPlot.pdf')
	else:
		pass

def Average(directory_to_search , plot_bool , show_plot):
	current = os.getcwd()
	os.chdir(directory_to_search)
	TheL = list()
	TheH = list()
	TheS = list()
	#Identify the Phi and Psi angels for each secondary structure
	with warnings.catch_warnings(record=True) as w:							#Supress Bio.PDB user warnings
		for TheFile in os.listdir('.'):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure(TheFile.split('.')[0] , TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				SS = list()
				Tor = list()
				for res in dssp:
					if res[2] == '-' or res[2] == 'T' or res[2] == 'S':		#Loop (DSSP code is - or T or S)
						SS.append('L')						#Secondary Structure
						Tor.append((res[4] , res[5]))				#(PHI , PSI)
					elif res[2] == 'G' or res[2] == 'H' or res[2] == 'I':		#Helix (DSSP code is G or H or I)
						SS.append('H')						#Secondary Structure
						Tor.append((res[4] , res[5]))				#(PHI , PSI)
					elif res[2] == 'B' or res[2] == 'E':				#Strand (DSSP code is B or E)
						SS.append('S')						#Secondary Structure
						Tor.append((res[4] , res[5]))				#(PHI , PSI)
				Tor = [(0.0 , x[1]) if x[0] == 360.0 else x for x in Tor]
				Tor = [(x[0] , 0.0) if x[1] == 360.0 else x for x in Tor]
				SecStruct = ''.join(SS)
				#Loop
				loop = list()
				for numb in SecStruct.replace('H' , '.').replace('S' , '.').split('.'):
					value = len(numb)
					if value == 0:
						pass
					else:
						loop.append(value)
				#PHI loop angels
				LPHI = [0]
				angel = 0
				for res , ang in zip(SS , Tor):			#Loop both lists

					if (res == 'H' or res == 'S'):		#If there is S or H from SS list ---> angel will become 0
						if angel != 0:			#And if the angel is a number other than 0
							LPHI.append(angel)	#Append it
							angel = 0		#Make the angel 0
					else:					#If there is L from SS list
						angel = angel + ang[0]		#Make angel value from Tor[0] list
				LPHI.append(Tor[-1][0])
				#PSI loop angels
				LPSI = list()
				angel = 0
				for res , ang in zip(SS , Tor):
					if (res == 'H' or res == 'S'):
						if angel != 0:
							LPSI.append(angel)
							angel = 0
					else:
						angel = angel + ang[1]
				LPSI.append(0)
				#Average loop angels
				LOOP = list()
				del LPHI[0]
				del LPSI[-1]
				for L , Ph , Ps in zip(loop , LPHI , LPSI):
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
				HPHI = list()
				angel = 0
				for res , ang in zip(SS , Tor):
					if (res == 'L' or res == 'S'):
						if angel != 0:
							HPHI.append(angel)
							angel = 0
					else:
						angel = angel + ang[0]
				#PSI helix angels
				HPSI = list()
				angel = 0
				for res , ang in zip(SS , Tor):
					if (res == 'L' or res == 'S'):
						if angel != 0:
							HPSI.append(angel)
							angel = 0
					else:
						angel = angel + ang[1]

				#Average helix angels
				HELIX = list()
				for L , Ph , Ps in zip(helix , HPHI , HPSI):
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
				SPHI = list()
				angel = 0
				for res , ang in zip(SS , Tor):
					if (res == 'L' or res == 'H'):
						if angel != 0:
							SPHI.append(angel)
							angel = 0
					else:
						angel = angel + ang[0]
				#PSI strand angels
				SPSI = list()
				angel = 0
				for res , ang in zip(SS , Tor):
					if (res == 'L' or res == 'H'):
						if angel != 0:
							SPSI.append(angel)
							angel = 0
					else:
						angel = angel + ang[1]
				#Average strand angels
				STRAND = list()
				for L , Ph , Ps in zip(strand , SPHI , SPSI):
					STRAND.append(((L) , (Ph/L) , (Ps/L)))
				#Organise final lists
				for item in LOOP:
					TheL.append(item)
				for item in HELIX:
					TheH.append(item)
				for item in STRAND:
					TheS.append(item)

			except Exception as TheError:
				print(TheFile , '<---' , TheError)
			print('Got average phi and psi angels and secondary structure lengths for\t' , TheFile.split('.')[0])
	os.chdir(current)
	TheLLe = list()
	TheLph = list()
	TheLps = list()
	TheHLe = list()
	TheHph = list()
	TheHps = list()
	TheSLe = list()
	TheSph = list()
	TheSps = list()
	for item in TheL:
		TheLLe.append(item[0])
		TheLph.append(item[1])
		TheLps.append(item[2])
	for item in TheH:
		TheHLe.append(item[0])
		TheHph.append(item[1])
		TheHps.append(item[2])
	for item in TheS:
		TheSLe.append(item[0])
		TheSph.append(item[1])
		TheSps.append(item[2])
	#Plot graph
	if plot_bool == 1:
		matplotlib.rcParams['axes.facecolor'] = '0.5'
		AX = mpl_toolkits.mplot3d.Axes3D(matplotlib.pyplot.figure())
		#L_plt = AX.scatter(TheLph , TheLps , TheLLe , c = '#038125' , s = 0.5)
		H_plt = AX.scatter(TheHph , TheHps , TheHLe , c = '#c82100' , s = 0.5)
		S_plt = AX.scatter(TheSph , TheSps , TheSLe , c = '#ffe033' , s = 0.5)
		AX.set_title('Number of Helices and Strands Plot' , y = 1.02)
		AX.set_xlim3d(-180, 180)
		AX.set_ylim3d(-180, 180)
		AX.set_xlabel('Phi')
		AX.set_ylabel('Psi')
		AX.set_zlabel('Length of Secondary Structure')
		#AX.legend([L_plt , H_plt , S_plt], ['Loop' , 'Helix' , 'Strand'] , markerscale = 10)
		AX.legend([H_plt , S_plt], ['Helix' , 'Strand'] , markerscale = 10)
		if show_plot == 1:
			matplotlib.pyplot.savefig('SSlenPlot.pdf')
			matplotlib.pyplot.show()
		else:
			matplotlib.pyplot.savefig('SSlenPlot.pdf')
	else:
 		pass

def KM(TheFile):
	''' Calculates the cluster centers of the Ramachandran plot, [output of the RamaPlot() function] '''
	''' Returns cluster centers '''
	df = pandas.read_csv(TheFile , sep = ';')
	L = df[['Loop PHI' , 'Loop PSI']]
	H = df[['Helix PHI' , 'Helix PSI']]
	S = df[['Strand PHI' , 'Strand PSI']]
	MLL = sklearn.cluster.KMeans(n_clusters = 5 , n_init = 10 , init = 'k-means++' , max_iter = 300).fit(L)
	MLH = sklearn.cluster.KMeans(n_clusters = 3 , n_init = 10 , init = 'k-means++' , max_iter = 300).fit(H)
	MLS = sklearn.cluster.KMeans(n_clusters = 5 , n_init = 10 , init = 'k-means++' , max_iter = 300).fit(S)
	#Loop plot
	matplotlib.rcParams['axes.facecolor'] = '0.8'
	matplotlib.pyplot.scatter(df['Loop PHI'] , df['Loop PSI'] , color = '#038125' , s = 0.5)
	matplotlib.pyplot.scatter(MLL.cluster_centers_[:,0] , MLL.cluster_centers_[:,1] , color = 'black')
	matplotlib.pyplot.title('Ramachandran Plot')
	matplotlib.pyplot.xlabel('Phi Angels')
	matplotlib.pyplot.ylabel('Psi Angels')
	matplotlib.pyplot.ylim(-180 , 180)
	matplotlib.pyplot.xlim(-180 , 180)
	matplotlib.pyplot.show()
	#Helix plot
	matplotlib.pyplot.scatter(df['Helix PHI'] , df['Helix PSI'] , color = '#c82100' , s = 0.5)
	matplotlib.pyplot.scatter(MLH.cluster_centers_[:,0] , MLH.cluster_centers_[:,1] , color = 'black')
	matplotlib.pyplot.title('Ramachandran Plot')
	matplotlib.pyplot.xlabel('Phi Angels')
	matplotlib.pyplot.ylabel('Psi Angels')
	matplotlib.pyplot.ylim(-180 , 180)
	matplotlib.pyplot.xlim(-180 , 180)
	matplotlib.pyplot.show()
	#Strand plot
	matplotlib.pyplot.scatter(df['Strand PHI'] , df['Strand PSI'] , color = '#ffe033' , s = 0.5)
	matplotlib.pyplot.scatter(MLS.cluster_centers_[:,0] , MLS.cluster_centers_[:,1] , color = 'black')
	matplotlib.pyplot.title('Ramachandran Plot')
	matplotlib.pyplot.xlabel('Phi Angels')
	matplotlib.pyplot.ylabel('Psi Angels')
	matplotlib.pyplot.ylim(-180 , 180)
	matplotlib.pyplot.xlim(-180 , 180)
	matplotlib.pyplot.show()
	#Print values
	print('Loop cluster centers:\n' , MLL.cluster_centers_ , '\n')
	print('Helix cluster centers:\n' , MLH.cluster_centers_ , '\n')
	print('Strand cluster centers:\n' , MLS.cluster_centers_ , '\n')
#--------------------------------------------------------------------------------------------------------------------------------------
#Database()
#RamaPlot('PDBDatabase' , 1)
#TheList = Numbers('PDBDatabase' , 1 , 1)
#Probability(TheList , 1 , 1)
#Length('PDBDatabase' , 1 , 1)
#Average('PDBDatabase' , 1 , 1)
KM('RamaPlot.csv')
