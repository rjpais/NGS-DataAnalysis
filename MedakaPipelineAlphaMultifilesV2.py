#============================================================================================================================
"""
MEDAKA PIPELINE VERSION ALPHA VERSION 2.1 Multifile analysis 
AUTHOR: RICARDO JORGE PAIS @ INSA   
DATE OF LAST UPDATE: 22/12/2020
THE SCRIPT PERFORMS AN AUTOMATED RUN OF MEDAKA METHOD FOR NGS SEQUENCING DATA USING OXFORD NANOPORE TECKNOLOGIES 
AND MANIPULATES THE OUTPUTS FOR INTEGRATION IN THE INSAFLU PLATFORM.  

AS INPUTS, THE SCRIPT TAKES:
	* A PATH FOR THE READS FILE IN FASTQ OR FASTQ.GZ (STRING)  
	* A PATH FOR THE REFERENCE GENOME SEQUENCE IN FASTA OR FASTA.GZ (STRING)
	* THE NAME OF THE MEDAKA MODEL TO BE USED

AS OUTPUTS, THE SCRIPT GENERATES THE FOLOWING FILES ON A NEW GENERATED FOLDER:
	*  Consensus fasta file with sample ID
	*  Bam file
	*  zipped Coverage file .depth.gz 
	*  VCF Variant file with coverages
"""
#============================================================================================================================
import os 
import gzip
import time
from sys import exit  
import shutil
from tkinter import *
from tkinter import filedialog # for easy runing on files and for lab users without using insaflu
from tkinter.filedialog import askdirectory 

#1 function for getting the sample ID name from filepath defined by user 
#--------------------------------------------------------------------------------------------------------------------------
def Get_Sample_IDname (filepath):
	""" filepath is the location of sample reads file
	returns the name of the file  
	"""
	name = filepath.split(".")[0].split("/")[-1]
	return name
#-------------------------------------------------------------------------------------------------------------------------- 


#2 function for generating a consensus fasta file  
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


#3 function for generating sample reads coverage   
#-------------------------------------------------------------------------------------------------------------------------- 
def CoverageExtraction(bam): 
	""" passes the commands of samtools for generating a .depth file 
	Requires the entery of the consensus probabilities file 
	Also condenses the file using gzip method """
	Output_file = bam.split("calls_to_draft")[0] + "reads_coverage.depth"
	commands =  "samtools depth -aa -d0 " + bam + " > " + Output_file 
	os.system(commands)
	exist_status = os.system(commands)
	if (exist_status != 0):
		print('Fail to run samtools commands\n please ensure that the tool is installed and run again the pipeline')
		exit(0)		
	print ("reads_coverage.depth was created")
	with open(Output_file, "rb") as initial:
		with gzip.open(Output_file + ".gz", "wb" ) as zipped:  
			zipped.writelines(initial)
			print ("reads_coverage.depth.gz was created")
	return Output_file
#-------------------------------------------------------------------------------------------------------------------------- 
	

#4 function for variant calling based on medaka consensus probabilities file  
#-------------------------------------------------------------------------------------------------------------------------- 
def VariantCalling_Medaka(probs, ref): 
	""" passes the commands of samtools for generating a .depth file 
	Requires the entery of the consensus probabilities file 
	Also """
	Output_file = probs.split("consensus_probs")[0] + "medaka_variant.vcf"
	commands =  "medaka variant --verbose " + ref + " " + probs + " " + Output_file
	os.system(commands)
	exist_status = os.system(commands)
	if (exist_status != 0):
		print('Fail to run medaka variant call commands\n please ensure that the tool is installed and run again the pipeline')
		exit(0)	
	print ("medaka_variant.vcf file was created") 
	return Output_file
#-------------------------------------------------------------------------------------------------------------------------- 
 


#5 function for adding information on a fasta file  
#-------------------------------------------------------------------------------------------------------------------------- 
def Add_SampleIDinfo_fasta(fastafile, info):
	""" fastafile is the fasta file path 
	info is a string with the information to add on the header 
	"""
	newfasta = ""
	fasta = open (fastafile, "r")
	for line in fasta:
		if line[0]== ">":
			newfasta = newfasta + line[0:-1] + " " + info + "\n"
		else:
			newfasta = newfasta + line	
	fasta.close()
	New = open(fastafile, "w" )
	New.write(newfasta)
	print(info , " was added to consensus header") 
	New.close()
#-------------------------------------------------------------------------------------------------------------------------- 
 

#6 function for getting the mutational positions and relevant information from VCF file  
#-------------------------------------------------------------------------------------------------------------------------- 
def Get_Variant_PositionsMutationScores(VCFpath):
	""" 
	VCFpath is the variant.vcf file path
	returns a list with all positions
	"""
	POSITIONS, MUTATIONS, SCORES = [ ], [ ],[ ] 
	vcf_file = open( VCFpath, "r" )
	for line in vcf_file:
		if line[0] !="#":
			POSITIONS.append(line.split("\t")[1])
			MUTATIONS.append(line.split("\t")[3] + "-->" + line.split("\t")[4])
			SCORES.append(line.split("\t")[7].split("pred_q=")[1].split(";")[0] )
	vcf_file.close() 
	return [POSITIONS, MUTATIONS, SCORES]
#-------------------------------------------------------------------------------------------------------------------------- 


#7 function for getting the coverage of each mutation identified  
#-------------------------------------------------------------------------------------------------------------------------- 
def Get_Variant_CoverageValues(positions, covfilepath): 
	"""
	positions is the list of all positions of mutations
	covfilepath is the path of the file containing the coverage of all reads
	"""
	Var_coverage = []
	coverage_file = open( covfilepath, "r" )
	for line in coverage_file:
		data = line.split("\t")
		if data[1] in positions:
			Var_coverage.append(data[2])
	coverage_file.close() 
	return Var_coverage
#-------------------------------------------------------------------------------------------------------------------------- 


#8 function for adding coverage values to variant file (VCF)  
#-------------------------------------------------------------------------------------------------------------------------- 
def Write_VCF_with_CoverageValues(coveragevalues, VCFpath): 
	"""
	coveragevalues is the list of coverage values for each mutation 
	VCFpath is the path of the file with variants
 	a new VCF file is produced that overwrites the previous one   
        """
	VCF2, n, k  = "", 0, 0
	addInfo = '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth value of the nucleotide position">\n'	 
	VCF_file = open( VCFpath, "r" )
	for line in VCF_file:
		if line[0] == "#" and k != 10:
			VCF2 = VCF2 + line
		elif line[0] == "#" and k == 10:
			VCF2 = VCF2 + addInfo 
		elif line[0] != "#": 
			C = line.split("\t")
			INFO_new = C[7] + ";DP=" + coveragevalues[n].split("\n")[0] 
			row = C[0]+"\t"+C[1]+"\t"+C[2] +"\t"+C[3]+"\t"+C[4]+"\t"+C[5]+"\t"+C[6]+"\t"+ INFO_new +"\t"+C[8]+"\t"+C[9] 
			VCF2 = VCF2 + row 
			n=n+1
		k = k+1		
	VCF_file.close()
	VCF_file = open(VCFpath, "w" )
	VCF_file.write(VCF2)
	VCF_file.close()
	print("medaka_variant.vfc was modified to have coverage values as DP=value")
#-------------------------------------------------------------------------------------------------------------------------- 


#9 function for renaming output files insaflu with insaflu standard   
#-------------------------------------------------------------------------------------------------------------------------- 
def INSaFlu_files_renamming(samplepath):
	""" 
	samplepath is a the reads file path
	"""
	path =   samplepath.split(".")[0]
	files = os.listdir(path)
	Prefix = path.split("/")[-1]
	for File in files:
		OldFilePath = path+"/"+File
		NewFile = Prefix+"."+File
		os.rename(OldFilePath, path+"/"+ NewFile)
		print ("file ", File , " was renamed to :", NewFile) 
#-------------------------------------------------------------------------------------------------------------------------- 


#10 function for removing files that are not required for storage   
#-------------------------------------------------------------------------------------------------------------------------- 
def UnecessaryFiles_remove(Gpath, Spath):
	""" 
	Gpath is the genome reference file path
	Spath is the sample reads reference file path
	"""
	Extensions = [".mmi", ".fai" ]
	for extension in Extensions:
		os.remove(Gpath+extension)
		print("file *", extension, " removed" )
	path =   Spath.split(".")[0]
	files = os.listdir(path)
	for File in files:
		if File.split(".")[-1] == "depth" or File.split(".")[-1] == "hdf":
			os.remove(path+"/"+File) 
			print("file ", File, " removed" )
#-------------------------------------------------------------------------------------------------------------------------- 

#10 function for filtering and get only high quality reads  
#-------------------------------------------------------------------------------------------------------------------------- 
def HQfilterReads(path, Q ):
	""" passes the commands of NanoFilt for generating a fastq file only with better quality reads  
	Requires the entery of full path of fastq reads file and the minimum quality to filter (Q)  
	Also """
	Output_file = path.split(".")[0] + "_HQonly.fastq.gz"	
	commands =  "gunzip -c " + path + " | NanoFilt -q " + str(Q) + "  | gzip > " + Output_file
	print ("\n ...filtering reads with quality > Q", str(Q), " \n ")
	os.system(commands)
	exist_status = os.system(commands)
	if (exist_status != 0):
		print('Fail to run NanoFilt tool commands for HQ reads filtering \n please ensure that the tool is installed and run again the pipeline')
		exit(0)		 
	return Output_file
	

# NANOPORE MEDAKA PIPELINE   
start = time.time()
print ("=========================================================================")
print ("AUTOMATED PIPELINE FOR MiniON READS MULTIPLE FILES PROCESSING   ")
print ("=========================================================================")


# Inputs that need to be passed by user on commandline usage

model = "r941_min_fast_g303"   
Tk().withdraw()
path = askdirectory(title = "open folder with  multiple sample reads files from MiniON" )  
Tk().withdraw()
RefGenome_path = filedialog.askopenfilename( title = "open reference genome fasta file" ) 
FILES = os.listdir (path)
Q = input("insert quality threshold filter value and press enter :   ") 

# AUTOMATED AND SYSTEMATIC DATA PROCESSING
N = 0 
for FileName in FILES:
	if FileName.split(".")[1] == "fastq":   
		N = N+1
		print("\n\n\n ...processing sample ", FileName.split(".")[0])
		sample_reads_path = path + "/" + FileName     
		HQsample_reads_path= HQfilterReads( sample_reads_path, int(float(Q)) )
		MedakaOutputs = Medaka_consensus_prediction (HQsample_reads_path , RefGenome_path , model)
		Consensus = MedakaOutputs[2]
		ProbFile = MedakaOutputs[1]
		BAMfile = MedakaOutputs[0]
		os.remove(path + "/" + FileName.split(".")[0] + "_HQonly.fastq.gz" )
		sampleIDname = Get_Sample_IDname (sample_reads_path)  
		Add_SampleIDinfo_fasta(Consensus, sampleIDname )     #  Manipulation of Consensus file header 
		SampleCoverageFile = CoverageExtraction(BAMfile)
		VCFfile = VariantCalling_Medaka(ProbFile, RefGenome_path)
		MutINFO = Get_Variant_PositionsMutationScores(VCFfile)
		Coverages = Get_Variant_CoverageValues(MutINFO[0], SampleCoverageFile)
		print ( "\nTotal of ", len(MutINFO[0]), " putative mutations detected!!")
		Write_VCF_with_CoverageValues(Coverages, VCFfile )
		UnecessaryFiles_remove(RefGenome_path, sample_reads_path)
		INSaFlu_files_renamming(sample_reads_path) 
		pTime = time.time() - start
		print ("\n=========================================================================")
print ("*******************END OF PROCESS*****THANK YOU *****************************   ")
print ("              Total number of samples processed = ", N )
print ("              Total processing time  = ", round(pTime , 1 ), " seconds ")
print ("              Processing time per sample  = ", round(pTime/N , 1 ), " seconds ")
print ("=========================================================================")

