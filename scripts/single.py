import os 

allo = open("all_annmax.csv",'r')

reps = open("all_annmax_single.csv", "w")


read = allo.readline().strip()

reps.write("%s\n" % read) 
	


while True:

	read = allo.readline().strip()

	if read == "":
		break

	row = read.split(";")

	

	if len(row[4])==1 and len(row[5])==1:
		reps.write("%s\n" % read)