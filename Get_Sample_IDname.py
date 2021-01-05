
#Function for getting the sample ID name from filepath defined by user 
#--------------------------------------------------------------------------------------------------------------------------
def Get_Sample_IDname (filepath):
	""" filepath is the location of sample reads file
	returns the name of the file  
	"""
	name = filepath.split(".")[0].split("/")[-1]
	return name
#-------------------------------------------------------------------------------------------------------------------------- 
