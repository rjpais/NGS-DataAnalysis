
# function for removing files that are not required for storage   
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
