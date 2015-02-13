import os
import glob

fPath = os.path.abspath(".")+"/listofeventnumbers"
outfile = open(fPath,"w")

for i in range(1,174):

	filePath = os.path.abspath(".")+"/Set_"+str(i)+"/out"
	infile = open(filePath,"r")
	input = infile.readlines()
	infile.close()

	for j in range(len(input)):
		if("Eventnumber:" in input[j]):
			number = input[j].split("Eventnumber: ")[1]
			outfile.write(number)
		else:
			continue

outfile.close()
