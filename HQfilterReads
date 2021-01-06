
#11 function for filtering and get only high quality reads  
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
#-------------------------------------------------------------------------------------------------------------------------- 
