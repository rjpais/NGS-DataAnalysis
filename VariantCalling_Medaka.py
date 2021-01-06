
# function for variant calling based on medaka consensus probabilities file  
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
 
