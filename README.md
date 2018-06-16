# AIDatabaseAnalyse
A script that uses Machine Learning to analyse all protein structures in the Protein DataBank (PDB) database

## Requirements:
1. Use the following command (in GNU/Linux) to install all necessary programs and python modules for this script to run successfully:

`sudo apt update && sudo apt install python3-pip python3-tk DSSP gnuplot && sudo pip3 install biopython numpy scipy pandas scikit-learn matplotlib`

## Description:
This is a script that uses Machine Learning to analyse all protein structures in the Protein DataBank (PDB) database. It primary focus is as follows:
* Find the ideal Helix/Strand/Loop torsion angels (average out the Ramachandran plot)
* Find the average length of Helixes and Stands
* Find the ratio of Helix to Stand

## How To Use:
1. Use the following command to run the script:

`python3 AIDatabaseAnalyse.py`

2. Computation time is very long (around 72 hours), and around +100GB of disk space is needed (mainly used to download and manipulate the PDB database).

**This script got replaced by the more advanced [ProtAI](https://github.com/sarisabban/AIDeNovo) script**
