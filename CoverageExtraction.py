
#Function for generating sample reads coverage   
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
	
