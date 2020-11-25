# 3DWSPA
3D chromosome reconstruction using Johnson's shortest path algorithm

# Contents of the folder

**src**: 3DWSPA.py is the source code<br/>
**parameters** : input_parameters.txt is the input parameter file to run the code<br/>
**examples**: contains example data and outputs generated from 3DWSPA for simulated Datasets<br/>

# Hi-C Input Data:

This code uses simulated data downloaded from https://umass.app.box.com/s/kbf2kjbk1sdpf4blceoajyo5lfyhx11s [1]

# Input matrix file format:

 This model acces N*N space delimited **Square Matrix** <br/>

# Usage

To run the code you execute the command : python3 3DWSPA.py input_parameters.txt <br/>

Input Parameters Contains <br/>
INPUT_FILE: Hi-C N*N contact matrix <br/>
OUTPUT_FOLDER: Output folder where pdb file and log files needs to be created.<br/>

# Output

3DWSPA produces 2 files <br/>

*.pdb: contains the model which can be visualized using pyMol <br/>
*_log.txt: Log files prints best average root means square error (RMSE), average correlation of Spearman's (SCC) and Pearson's correlations(PSC) of separate models generated using different converversion factor , which varies from 0 to 2.

# Citation
[1] Zou, Chenchen, Yuping Zhang, and Zhengqing Ouyang. "HSA: integrating multi-track Hi-C data for genome-scale reconstruction of 3D chromatin structure." Genome biology 17, no. 1 (2016): 40.
