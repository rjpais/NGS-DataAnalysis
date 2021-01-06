#function for getting the coverage of each mutation identified  
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
