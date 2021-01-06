
# function for adding coverage values to variant file (VCF)  
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
