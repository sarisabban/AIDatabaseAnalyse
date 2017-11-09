#!/usr/bin/python3

import os , math , gzip , warnings , random , Bio.PDB , matplotlib.pyplot , mpl_toolkits.mplot3d
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






def DBSCAN(TheList):
	print(TheList)















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
#TheList = ([0.0, -74.0, -75.0, -142.0, -74.7, -87.3, 69.4, -97.9, -103.1, -114.6, -85.7, -66.5, 100.5, -81.3, -44.5, -108.6, 68.2, -142.7, -113.4, -67.3, 70.5, -96.9, 76.8, 77.3, -118.3, -132.4, 83.6, -86.2, -153.3, -155.7, -89.5, -115.8, -83.2, -113.6, -109.0, -63.6, -140.7, 72.9, -114.2, -134.3, -65.9, -130.0, -127.9, -51.9, 84.6, -60.8, -89.4, -58.3, -73.5, -123.3, -88.6, -59.5, -70.7, -84.1, -103.4, -70.4, -64.8, -118.1, -72.2, -110.1, 72.9, -97.6, -165.4, 0.0, 0.0, -64.9, -90.1, 62.0, -143.4, -79.0, -124.8, -61.7, -124.9, -86.7, 98.4, -99.3, 80.9, -82.8, -55.9, -69.9, -96.6, 39.1, -87.3, 102.6, -92.5, -66.6, -91.9, -98.4, -57.3, -103.8, -78.6, -107.7, -88.4, -65.2, -80.5, -132.1, -64.7, -30.8, -85.9, -130.3, -33.9, -80.8, -171.6, -96.1, -66.1, -89.7, -115.0, 52.9, -96.1, -113.0, -57.7, 101.2, -76.1, -120.2, -49.5, -77.2, -47.0, -59.1, -139.7, -110.5, -88.5, -55.4, -119.8, -122.3, -107.2, -74.9, -97.4, -74.4, -90.0, 65.4, -117.7, -85.2, -164.8, 89.4, -74.6, -66.8, -73.1, -112.9, -58.7, -56.1, -66.8, -148.0, -88.3, -119.3, -72.0, -89.6, -94.1, -79.8, 87.3, 64.9, -121.4, 121.5, -70.6, -72.5, -103.3, -58.7, -84.5, 99.4, -127.2, -52.2, -55.2, -145.0, -44.4, -147.8, -86.2, -72.1, -105.0, -148.0, -91.8, -121.6, -139.8, -80.5, -102.2, -40.6, -68.7, 85.9, -90.4, 131.3, -114.1, -54.0, -56.5, -82.8, 59.4, -105.4, -71.0, -97.2, 90.4, -137.1, -79.1, -59.9, 75.9, -100.8, -93.7, -60.2, -74.8, -96.0, 43.2, -104.0, -69.7, -135.5, -108.5, -137.9, -80.3, -107.0, -119.5, -128.8, -56.9, -114.3, -135.4, -109.9, -56.6, 71.8, -67.8, 137.0, 0.0, 0.0, -82.9, -90.3, -111.1, -58.7, -99.3, -73.7, -136.2, -62.3, 96.7, -67.7, -98.9, -76.7, -160.0, -44.6, -124.3, 57.4, -114.5, -88.0, -115.9, -56.9, 64.8, -126.4, -64.0, -56.5, 54.4, 68.1, -127.6, -76.6, -55.9, 83.2, -95.0, -52.2, -67.5, -75.1, 75.6, -121.1, -93.0, 51.9, -51.0, -112.2, -68.0, -177.1, -134.0, -86.8, -140.8, -88.0, -76.7, -63.8, -78.6, -81.0, 94.3, -83.3, -144.7, -70.3, -62.0, -80.2, -71.7, -65.5, -53.0, -80.8, -69.2, -112.6, 80.2, -73.5, -60.6, -53.5, -159.8, 50.1, 77.0, -130.7, -87.5, 89.3, -87.0, -73.5, -52.5, -53.8, -83.2, -98.7, -119.8, 54.9, -95.0, -70.8, -92.2, -131.0, -113.5, -62.8, -75.0, -102.7], [-132.8, 162.2, 57.1, 165.1, -56.5, 4.6, 31.9, 148.6, 143.4, 121.6, 163.8, -28.3, -37.6, 130.6, -50.7, 13.4, 26.0, 164.9, 163.2, 157.7, -127.4, 5.4, 111.2, 39.7, 46.9, 35.6, 165.6, -18.8, -176.8, 68.0, -4.7, 96.6, -22.2, -115.4, 177.4, -32.7, 123.6, -161.0, 15.8, 145.1, 145.5, 172.3, 85.0, 128.0, 4.7, -22.0, -1.6, -41.0, -4.7, 86.7, -10.6, 128.3, 156.7, 8.8, 133.6, 151.6, -50.5, 162.7, -9.6, 13.2, 4.4, 150.3, 0.0, 0.0, -62.1, -59.3, -2.9, 45.5, 143.0, 140.8, 148.0, 159.5, -0.3, 125.5, 14.4, -25.7, -1.8, -37.9, 128.6, 124.3, 7.3, -106.5, -25.5, -1.6, -45.0, 134.8, 146.2, 156.7, 127.5, 78.3, 127.1, 119.0, 115.8, -28.9, -7.1, 74.3, -51.3, -55.6, 162.5, 149.6, 129.6, -37.2, 158.8, 49.3, -52.2, -23.3, 11.0, 86.5, 165.1, 179.2, 122.5, -31.9, 133.8, 77.0, -36.4, -7.6, 112.8, -38.7, 145.6, 98.6, -174.6, -13.9, -8.2, 119.3, 156.8, 1.7, 153.0, -163.0, 134.6, 36.9, -0.4, 89.9, -177.5, 2.0, 132.6, 116.2, 117.9, -26.1, 136.1, 153.2, -30.8, 152.4, 146.7, 71.2, -14.3, -43.5, 159.9, 179.4, -59.5, 77.9, 123.1, -36.9, 147.6, 144.5, -166.4, -37.1, -12.4, 16.5, 117.4, 145.5, -56.7, 149.2, 146.8, 166.5, 141.0, -24.0, 23.2, -61.9, 112.5, 18.9, -14.7, 113.0, -179.1, -36.9, 1.3, -42.8, 154.8, 136.8, 12.1, -40.7, -56.3, -10.8, 53.0, -49.8, 142.9, 168.7, 23.4, 81.4, 153.6, -11.0, 45.7, -31.6, 147.0, 169.6, 129.3, 15.5, 52.2, 145.5, -29.2, 146.2, -13.7, 134.5, 88.3, -169.8, 128.2, 153.9, 160.3, 26.0, -25.6, 107.0, 153.3, 27.2, 105.3, 0.0, 0.0, 135.1, 126.5, 124.9, 94.6, -26.5, 110.5, -35.8, 152.1, 140.2, -21.2, 153.5, 122.4, -4.4, 147.0, 126.0, 18.3, -129.2, 24.0, 64.4, 158.4, 132.0, 19.7, 162.7, 139.2, 158.3, 32.8, -38.4, 1.0, 151.8, 138.4, 3.5, 127.3, 137.2, -5.9, -9.2, -98.5, -18.1, 141.5, 55.9, 136.2, 162.4, 112.2, -178.3, 35.3, -37.5, 150.9, 146.0, 148.3, 135.3, 167.5, -3.5, 150.4, 137.2, 178.4, 151.1, 144.6, 140.8, 151.8, 162.8, 141.0, 152.5, -8.9, 13.2, 4.8, 178.7, -37.7, 120.8, 154.2, 42.4, 9.2, -11.4, 144.4, -9.6, 158.5, 160.0, 123.2, 142.0, -15.2, -42.2, 7.4, 18.8, 175.4, -0.2, 2.3, 166.2, -59.3, -12.9, -3.7, 0.0], [34.8, -81.1, -111.1, -53.4, -63.2, -67.7, -71.9, -62.3, -63.0, -61.1, -59.7, -62.7, -61.7, -68.3, -55.7, -64.3, -69.3, -54.5, -63.9, -53.5, -63.5, -85.1, -51.2, -69.4, -68.0, -62.7, -59.5, -65.2, -64.6, -65.4, -64.4, -74.3, -75.3, -94.9, -61.7, -60.6, -66.0, -59.0, -66.2, -68.1, -62.0, -62.8, -60.5, -68.0, -60.3, -63.6, -67.1, -60.5, -58.7, -61.8, -68.5, -63.4, -63.1, -66.3, -62.2, -60.8, -61.5, -64.0, -55.4, -91.3, -99.1, -60.8, -54.8, -68.7, -85.0, -55.2, -68.2, -62.7, -60.3, -66.1, -67.4, -81.3, -48.8, -51.9, -66.6, -74.0, -85.1, -48.6, -52.0, -74.9, -66.5, -64.8, -56.7, -62.4, -62.5, 47.8, 61.7, 70.2, -64.1, -65.8, -61.4, -62.7, -69.3, -67.5, -56.4, -59.0, -56.7, -65.3, -59.7, -69.6, -80.9, -65.3, -72.9, -68.8, -52.6, -71.8, -65.1, -61.3, -65.2, -58.7, -61.1, -63.0, -64.1, -67.9, -63.1, -34.7, -63.2, -69.7, -62.3, -61.5, -63.7, -64.7, -94.0, -63.5, -52.7, -94.4, -54.4, -56.7, -44.0, -69.7, -73.6, -76.5, -89.5, -80.9, -61.1, -56.0, -73.4, -52.8, -52.8, -81.0, -78.1, -71.6, -48.2, -79.2, -62.6, -67.1, -90.4, -89.6, -62.8, -59.6, -70.5, -76.1, -87.2, -60.7, -95.9, -67.8, -73.0, -2.4, -57.3, -48.8, -70.4, -77.1, -75.8, -75.9, -69.8, -78.1, -64.9, -69.5, -73.7, -71.7, -62.9, -88.1, -66.1, -82.3, -64.2, -58.2, -73.8, -74.9, -80.7, -76.0, -57.6, -128.9, -72.6, -77.3, -76.8, -64.5, -52.2, -89.5, -57.5, -46.9, -62.4, -69.9, -69.2, -66.7, -72.0, -91.6, -111.1, -73.3, -54.1, -61.4, -104.5, -56.7, -67.5, -66.8, -60.1, -66.7, -48.6, -73.8, -63.5, -59.4, -58.9, -102.6, -91.5, -58.4, -74.5, -75.5, -90.1, -61.7, -60.7, -66.3, -46.6, -78.7, -63.8, -57.1, -87.0, -75.6, -65.8, -67.0, -69.2, -50.6, -46.9, -71.9, -72.2, -76.0, -92.5, -51.0, -63.5, -98.4, -65.7, -52.6, -66.6, -87.8, -55.3, -54.9, -65.9, -68.8, -57.4, -58.4, -66.3, -62.8, -58.5, -76.0, -132.2], [-126.0, 7.9, -9.1, -50.8, -45.5, -29.3, -42.8, -38.4, -50.7, -46.6, -52.4, -29.9, -29.3, -54.9, -40.9, -43.8, -11.6, -49.7, -48.2, -43.1, -30.8, -10.5, -44.6, -29.6, -40.3, -45.9, -46.9, -35.8, -47.7, -39.2, -35.1, -46.9, -30.9, -4.1, -38.6, -47.5, -37.4, -40.5, -39.8, -37.7, -43.5, -43.7, -39.3, -48.8, -41.5, -35.0, -44.1, -45.5, -40.7, -47.4, -33.3, -39.5, -39.2, -40.6, -48.6, -41.0, -41.1, -53.4, -48.8, -22.0, -49.5, -56.6, -54.6, -51.8, -3.1, -40.4, -48.4, -40.6, -45.7, -33.6, -38.2, -14.2, -48.9, -41.7, -33.2, -32.1, -27.0, -44.9, -38.1, -32.4, -39.1, -48.2, -44.1, -45.9, -25.1, 39.3, 25.2, 22.4, -26.0, -32.0, -48.1, -40.6, -40.4, -43.1, -48.1, -46.8, -44.1, -45.3, -44.6, -22.6, -27.6, -35.9, -40.6, -51.8, -40.0, -42.0, -41.2, -46.8, -39.4, -48.8, -47.4, -41.0, -49.1, -40.2, -16.7, -68.9, -36.7, -40.6, -44.4, -45.7, -44.7, -21.6, -13.7, -32.2, -40.0, -57.2, -41.0, -64.7, -46.6, -44.2, -18.8, -28.6, -9.7, -49.6, -45.5, -33.3, -52.3, -61.5, -36.1, -24.9, -39.8, -47.4, -38.9, -25.8, -43.9, -18.2, -29.2, -53.1, -41.9, -39.7, -31.5, -18.7, -39.5, -27.2, -24.7, -40.2, -66.9, -61.3, -58.2, -32.8, -51.4, -20.1, -35.8, -24.9, -31.0, -40.6, -39.8, -31.4, -29.8, -41.9, -20.2, -34.2, -33.7, -22.3, -53.2, -35.9, -39.8, -20.3, -30.9, -66.1, -27.6, -13.4, -19.5, -15.3, -32.7, -34.9, -37.0, -49.9, -52.8, -57.8, -26.5, -48.7, -32.6, -46.6, -40.5, 11.6, -21.6, -47.7, -51.0, -51.0, -17.6, -37.6, -46.7, -28.4, -59.1, -37.6, -56.1, -31.1, -39.7, -50.0, -28.2, -12.4, -6.8, -17.4, -11.4, -6.1, -18.0, -42.9, -46.6, -44.6, -54.5, -18.4, -50.1, -46.6, -47.2, -31.5, -42.7, -39.3, -45.3, -60.1, -61.6, -25.3, -21.6, -1.6, -13.4, -48.4, -25.5, -28.2, -27.2, -41.5, -10.4, -7.4, -38.9, -52.0, -33.0, -32.6, -43.4, -50.0, -36.7, -46.0, -34.0, 16.2, 34.7], [-111.6, -143.3, -111.8, -102.7, -161.2, -151.4, -104.4, -104.4, -106.8, -79.8, -124.2, -132.9, -133.5, -127.3, -84.1, -106.6, -122.4, -132.2, -98.8, -104.3, -85.9, -127.1, -95.3, -81.2, -140.4, -152.0, -73.2, -154.9, -119.1, -137.3, -129.4, -81.4, -62.5, -132.3, -136.5, -142.3, -132.9, -115.1, -161.5, -142.1, -73.4, -83.1, -124.4, -161.2, -155.1, -77.0, -108.7, -112.9, -133.2, -59.2, -125.7, -140.5, -109.2, -81.8, -84.7, -65.0, -58.0, -82.9, -116.2, 171.8, -119.0, -78.1, -130.6, -139.1, -77.8, -100.0, -164.1, -69.0, -135.5, -69.6, -138.3, -124.4, -135.8, 160.3, -63.2, -77.0, -137.2, -113.5, -146.0, 176.4, -122.5, -114.2, -61.9, -106.4, -99.2, -136.5, -147.1, -86.7, -141.5, -121.4, -126.9, -119.5, -100.5, -151.3, -116.5, -97.1, -136.0, -131.8, -120.2, -124.0, -111.0, -141.1, -137.5, -75.6, -112.9, -134.3, -141.9, -102.1, -77.9, -94.8, -140.8, -140.9, -155.1, 112.6, -173.1, -112.0, -126.9, -145.7, -102.7, -117.2, -108.3, -92.9, -128.8, -144.3, -92.8, -154.8, -105.9, -128.4, -121.6, -132.1, -99.8, -86.5, -121.2, -88.8, -143.6, -98.1, -100.4, -141.8, -119.6, -133.0, -117.2, -148.7, -109.0, -109.7, -112.1, -129.6, -107.1, -126.5, -100.7, 48.0, -120.0, -163.6, -132.2, -116.8, -133.0, -109.0, -135.0, -116.1, -114.7, -55.9, -135.9, -95.2, -123.8, -112.0, -135.1, -123.1, -123.6, -141.7, -152.9, -120.2, -133.5, -131.4, -98.1, -111.5, -129.4, -91.3, -92.6, -124.3, -133.9, -109.9, -119.9, -111.9, -111.6, -140.6, -166.1, -80.3, -97.2, -73.7, -85.0, -114.9, -127.9, -159.6, -154.1, -116.1], [148.2, 123.1, 122.8, 118.6, 169.9, 158.9, 136.4, 120.7, 98.8, 129.0, 139.6, 148.2, 130.0, 94.4, 118.8, 125.7, 155.5, 169.3, 141.0, 156.4, 147.0, 121.5, 140.6, 145.5, 168.7, 161.0, 118.0, -178.8, 161.5, 157.1, 149.0, 80.8, 143.6, 124.2, 175.7, 120.4, 115.6, 163.3, 133.2, 148.0, 85.0, 108.3, 165.1, 155.1, 135.8, 117.6, 156.2, 115.9, 110.4, 155.1, 133.7, 132.9, 102.7, 113.2, 154.4, 141.9, 134.6, 131.9, -175.2, 169.6, 129.7, 113.9, 154.2, 162.4, 123.9, 96.4, 166.1, 144.2, 108.9, 128.1, 152.2, 112.2, 178.4, 167.9, 129.8, 163.2, 156.0, 125.9, -178.4, 178.1, 136.5, 130.2, 133.9, 108.4, 118.5, 152.0, 124.8, 138.1, 173.0, 141.7, 132.9, 130.2, 156.0, 124.6, 131.8, 145.0, 146.6, 155.2, 139.6, 147.6, 136.5, 135.1, 140.8, 125.7, -47.7, 145.6, 161.8, 130.5, 138.4, 128.9, 162.4, 156.7, 173.7, 169.7, 159.9, 119.5, 153.6, 132.7, 123.9, 127.2, 116.1, 142.8, 140.1, 144.1, 165.8, 139.4, 147.3, 154.5, 157.5, 124.4, 124.5, 125.6, 154.3, 136.9, 132.6, 112.4, 136.4, 135.9, 144.7, 129.8, 165.8, 106.5, 129.9, 117.7, 143.3, 129.7, 136.9, 90.8, 150.1, 63.1, 166.7, 138.3, 135.8, 137.9, 138.5, 144.5, 131.7, 119.4, 132.9, 139.4, 127.5, 137.4, 132.0, 140.2, 134.7, 142.1, 157.9, 158.0, 148.6, 144.0, 123.7, 131.0, 119.1, 121.7, 179.1, 161.1, 128.4, 130.0, 133.3, 117.5, 125.6, 126.0, 124.5, 154.6, 141.7, -16.3, 128.6, 152.4, 123.3, 128.1, 139.0, 175.2, 160.5, 131.0])
#DBSCAN(TheList)

Numbers('PDBDatabase' , 1 , 1)

'''
* DBSCAN center of clusters
* average angles plot
* ratio of helix to strand to loop to each other and to size of protein
* leangth of helix and strand to size of protein
'''
