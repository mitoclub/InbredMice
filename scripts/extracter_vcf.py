import os 



#files = ["Control-Polg-1","Control-Polg-2","Control-Polg-3","Parkin-Polg-1","Parkin-Polg-2","Parkin-Polg-3"]
#files = ["MGI-1","MGI-2","MGI-3","MGI-4","MGI-97","MGI-98","MGI-99","MGI-100","MGI-101","MGI-102","MGI-103","MGI-104"]
files = ["MitoLab"]

#prepath = "C:\\Users\\Abathur\\Desktop\\ParkPolGbam\\first_snp\\"
#prepath = "C:\\Users\\Abathur\\Desktop\\ParkPolGbam\\second_snp\\"
prepath = "D:\\________ESCAPE\\Science_122021\\InbredMice\\MtBAMs\\VCFs\\HaplotypeCaller\\"


#chrs=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","21","22","X","Y","MT"]
chrs=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr21","chr22","chrX","chrY","chrM"]


#reps = open(prepath + "first_snp.csv" , "w")
reps = open(prepath + "snps.csv" , "w")






files=os.listdir("HaplotypeCaller")





for fi in files:


	allo = open(prepath + "%s" % fi,'r')

	arr = []



	while True:

		read = allo.readline().strip()

		if "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" in read:
			break

	count = 0
	sumvaf = 0.0

	while True:

		read = allo.readline().strip()

		if read == "":
			break

		row = read.split("	")
		info = row[7].split(";")

		if row[0] in chrs:

			af = info[1].replace("AF=","").split(",")
			for a in af:
				count += 1
				sumvaf += float(a)

	if count != 0:
		reps.write("%s;%s;%s\n" % ( fi,count,sumvaf/count ))
	else:
		reps.write("%s;%s;%s\n" % ( fi,count,"none" ))