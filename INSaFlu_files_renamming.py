
# function for renaming output files insaflu with insaflu standard   
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
