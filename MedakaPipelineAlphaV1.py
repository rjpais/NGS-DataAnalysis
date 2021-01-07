#============================================================================================================================
"""
MEDAKA PIPELINE VERSION ALPHA 1
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
	I, M, R  = samplepath, model, refpath
	O = samplepath.split(".")[0] # output folder
	if M == "default":
 		commands =  "medaka_consensus -i "+ I +" -d "+ R +  " -o " + O
	else:
 		commands =  "medaka_consensus -i "+ I +" -d "+ R +  " -o " + O + " -m " + M 
	os.system(commands)
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
	print ("medaka_variant.vcf file was created") 
	os.system(commands)
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
 

#6 function for getting the mutational positions from consensus  
#-------------------------------------------------------------------------------------------------------------------------- 
def Get_Variant_Positions(VCFpath):
	""" 
	VCFpath is the variant.vcf file path
	returns a list with all positions
	"""
	mutpos = []
	vcf_file = open( VCFpath, "r" )
	for line in vcf_file:
		if line[0] !="#":
			position = line.split("\t")[1]
			mutpos.append(position)
	vcf_file.close() 
	return mutpos
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
	VCF2, n  = "", 0
	header2 = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tCOVERAGE\tFORMAT\tSAMPLE\n"	 
	VCF_file = open( VCFpath, "r" )
	for line in VCF_file:
		if line[0:2] == "##":
			VCF2 = VCF2 + line
		elif line[0] == "#" and line[1] != "#":
			VCF2 = VCF2 + header2
		else: 
			C = line.split("\t")
			Kn = "DP=" + coveragevalues[n].split("\n")[0] 
			row = C[0]+"\t"+C[1]+"\t"+C[2] +"\t"+C[3]+"\t"+C[4]+"\t"+C[5]+"\t"+C[6]+"\t"+C[7]+"\t"+ Kn +"\t"+C[8]+"\t"+C[9] 
			VCF2 = VCF2 + row 
			n=n+1		
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
	 

# NANOPORE MEDAKA PIPELINE   
start = time.time()
print ("=========================================================================")
print ("AUTOMATED RUN OF MEDAKA METHOD FOR A ONT READS (ADAPTED FOR INSaFlu)  ")
print ("=========================================================================")


# 0 Inputs that need to be passed by user on the platform insaflu
sample_reads_path = "/home/ricardo/TELE-VIRproject/NANOporeTESTdata/ERR4082025_1.fastq.gz"
RefGenome_path = "/home/ricardo/TELE-VIRproject/NANOporeTESTdata/SARS_CoV_2_Wuhan_Hu_1_MN908947.fasta"
model = "r941_min_fast_g303"   


# STEP 1 runing medaka for generating consensus prediction, BAM and probabilities files   
print("\n\n STEP 1: Running Medaka consensus prediction tool")
print ("=========================================================================")
MedakaOutputs = Medaka_consensus_prediction (sample_reads_path , 
RefGenome_path , model)
Consensus = MedakaOutputs[2]
ProbFile = MedakaOutputs[1]
BAMfile = MedakaOutputs[0]
sampleIDname = Get_Sample_IDname (sample_reads_path)  
Add_SampleIDinfo_fasta(Consensus, sampleIDname )     #  Manipulation of Consensus file header 
print("Key output files generated with medaka_consensus method")
for File in MedakaOutputs:
	print (File)
print ("\n=========================================================================")


# STEP2 Generation of a coverage/depth file from medaka BAM file
print("\n\n STEP 2: Computing reads coverage using samtools")
print ("=========================================================================\n")
SampleCoverageFile = CoverageExtraction(BAMfile)
print ("=========================================================================")


# STEP3 Variant calling  from medaka probabilities and reference genome with insaflu format
print("\n\n STEP 3: Computing variant calls using medaka variant method")
print ("=========================================================================\n")
VCFfile = VariantCalling_Medaka(ProbFile, RefGenome_path)
Mutationalpositions = Get_Variant_Positions(VCFfile)
Coverages = Get_Variant_CoverageValues(Mutationalpositions, SampleCoverageFile)
print("Position Coverage" )
for i, pos in enumerate(Mutationalpositions):
	print(pos, "\t", Coverages[i][0:-1] )	
print ( "Total of ", len(Mutationalpositions), " putative mutations detected!!")
#Write_VCF_with_CoverageValues(Coverages, VCFfile )


# STEP4 Removing and Renamming output files with insaflu standard form
print("\n\n STEP 4: Removing intermediate files and renaming output files")
print ("=========================================================================\n")
#UnecessaryFiles_remove(RefGenome_path, sample_reads_path)
INSaFlu_files_renamming(sample_reads_path) 
print ("=========================================================================")
print ("*******************END OF PROCESS*****THANK YOU *****************************   ")
print ("              Total processing time = ", round(time.time() - start , 1 ), " seconds ")
print ("=========================================================================")

