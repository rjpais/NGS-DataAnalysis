# Requirements: 
import os 
from sys import exit  
import shutil

#Function for generating a consensus fasta file  
#-------------------------------------------------------------------------------------------------------------------------- 
def Medaka_consensus_prediction(samplepath ,refpath, model):
	""" passes the comamnds for running medaka consensus method for consensus sequence prediction.
	Requires the entery of the sample reads path (samplepath) and refrence sequence (refpath)  + the indication of the model to be used (model)
	if model = "default" it runs as original default model 
	Saves otput files on a folder with sampleIDname and return the 3 key outputs file paths for further processing      	
	"""
	I, M, R  = samplepath , model, refpath
	O = samplepath.split("_HQonly")[0]  # output folder
	output_exists = os.path.isdir(O)
	if output_exists == True:
		shutil.rmtree(O)
	if M == "default":
 		commands =  "medaka_consensus -i "+ I +" -d "+ R +  " -o " + O
	else:
 		commands =  "medaka_consensus -i "+ I +" -d "+ R +  " -o " + O + " -m " + M 
	os.system(commands)
	exist_status = os.system(commands)
	if (exist_status != 0):
		print('Fail to run medaka consensus tool commands\n please ensure medaka is installed and run again the pipeline')
		exit(0)		
	bamFile = O + "/calls_to_draft.bam" 
	ProbsFile = O + "/consensus_probs.hdf"
	Consensus = O + "/consensus.fasta"
	return [bamFile, ProbsFile, Consensus ]
#-------------------------------------------------------------------------------------------------------------------------- 
