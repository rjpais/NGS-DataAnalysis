function for adding information on a fasta file  
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
 
