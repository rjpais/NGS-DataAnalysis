#function for getting the mutational positions and relevant information from VCF file  
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
