import os 

#files = ["Control-Polg-1","Control-Polg-2","Control-Polg-3","Parkin-Polg-1","Parkin-Polg-2","Parkin-Polg-3"]
#files = ["MGI-1","MGI-2","MGI-3","MGI-4","MGI-97","MGI-98","MGI-99","MGI-100","MGI-101","MGI-102","MGI-103","MGI-104"]
files=os.listdir("data\\HC\\")

#prepath = "C:\\Users\\Abathur\\Desktop\\ParkPolGbam\\mutect_first\\"
prepath = "C:\\Users\\Abathur\\Desktop\\Science\\InbredMice\\data\\HC\\"
prepath2 = "C:\\Users\\Abathur\\Desktop\\Science\\InbredMice\\data\\HC_split\\"

#chrs=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr21","chr22","chrX","chrY","chrM"]
chrs=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","21","22","X","Y","MT"]






def mostwanted(anns):
	maxid=0
	maxscore=0
	
	for i in range(len(anns)):
		score=0
		row = anns[i].split("|")

		if "WARNING_TRANSCRIPT_NO_START_CODON" in row[15]:
			score+=20
		if "WARNING_TRANSCRIPT_INCOMPLETE" in row[15]:
			score+=15
		if "WARNING_TRANSCRIPT_NO_STOP_CODON" in row[15]:
			score+=10			

		if "HIGH" in row[2]:
			score+=5
		if "MODERATE" in row[2]:
			score+=3
		if "LOW" in row[2]:
			score+=2
		if "MODIFIER" in row[2]:
			score+=1



		if score > maxscore:
			maxid = i

	return anns[maxid],maxscore






reps5 = open(prepath2 + "count.csv", "w")
reps5.write("Name;muts;vaf\n") 




reps6 = open(prepath2 + "all_annmax.csv", "w")
reps6.write("File;Name;Chr;Pos;Ref;Alt;SNP;alt_allele_frequency;AN;BaseQRankSum;DP;ExcessHet;FS;InbreedingCoeff;MLEAC;MLEAF;MQ;MQRankSum;QD;ReadPosRankSum;SOR;Allele;Annotation;Annotation_Impact;Gene_Name;Gene_ID;Feature_Type;Feature_ID;Transcript_BioType;Rank;HGVS.c;HGVS.p;cDNA.pos/cDNA.length;CDS.pos/CDS.length;AA.pos/AA.length;Distance;ERRORS/WARNINGS/INFO;Score;lof\n") 
	





for fi in files:


	# onlycommon


	allo = open(prepath + fi,'r')

	arr = []


	reps = open(prepath2 + "%s_ann.csv" % fi, "w")
	# Functional annotations:   'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
	reps.write("Chr;Pos;Ref;Alt;SNP;Allele;Annotation;Annotation_Impact;Gene_Name;Gene_ID;Feature_Type;Feature_ID;Transcript_BioType;Rank;HGVS.c;HGVS.p;cDNA.pos/cDNA.length;CDS.pos/CDS.length;AA.pos/AA.length;Distance;ERRORS/WARNINGS/INFO;lof\n") 
	
	reps2 = open(prepath2 + "%s_lof.csv" % fi, "w") 
	# Predicted loss of function effects for this variant. Format:  'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'
	reps2.write("Chr;Pos;Ref;Alt;Gene_Name;Gene_ID;Number_of_transcripts_in_gene;Percent_of_transcripts_affected\n") 
	
	reps3 = open(prepath2 + "%s_nmd.csv" % fi, "w") 
	# Predicted nonsense mediated decay effects for this variant. Format:  'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'
	reps3.write("Chr;Pos;Ref;Alt;Gene_Name;Gene_ID;Number_of_transcripts_in_gene;Percent_of_transcripts_affected\n") 



	reps4 = open(prepath2 + "%s_annmax.csv" % fi, "w")
	reps4.write("Chr;Pos;Ref;Alt;SNP;alt_allele_frequency;AN;BaseQRankSum;DP;ExcessHet;FS;InbreedingCoeff;MLEAC;MLEAF;MQ;MQRankSum;QD;ReadPosRankSum;SOR;Allele;Annotation;Annotation_Impact;Gene_Name;Gene_ID;Feature_Type;Feature_ID;Transcript_BioType;Rank;HGVS.c;HGVS.p;cDNA.pos/cDNA.length;CDS.pos/CDS.length;AA.pos/AA.length;Distance;ERRORS/WARNINGS/INFO;Score;lof\n") 
	


	name = "_".join(fi.split("_phased_")[0].split("_")[1:])

	##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
	##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">

	##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
	##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
	##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
	##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
	##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
	##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
	##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
	##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
	##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
	##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
	##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
	##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">








	while True:

		read = allo.readline().strip()

		if "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" in read:
			break


	muts = 0
	sumass = 0.0
	while True:

		read = allo.readline().strip()

		if read == "":
			break

		row = read.split("	")
		info = row[7].split(";")
		assb = info[1].replace("AF=","").split(",")
		allas = 0
		vaf = []


		# print(fi)
		# print(info) 




		an = ""
		bqrs = ""
		dp = ""
		fs = ""
		eh = ""
		ic = ""
		mleac = ""
		mleaf = ""
		mq = ""
		mqrs = ""
		qd = ""
		rprs = ""
		sor = ""



		




		for inf in info:
			if  "AN=" in inf:
				an = inf.replace("AN=","")
			elif  "BaseQRankSum=" in inf:
				bqrs = inf.replace("BaseQRankSum=","")
			elif  "DP=" in inf:
				dp = inf.replace("DP=","")
			elif  "FS=" in inf:
				fs = inf.replace("FS=","")
			elif  "ExcessHet=" in inf:
				eh = inf.replace("ExcessHet=","")
			elif  "InbreedingCoeff=" in inf:
				ic = inf.replace("InbreedingCoeff=","")
			elif  "MLEAC=" in inf:
				mleac = inf.replace("MLEAC=","")
			elif  "MLEAF=" in inf:
				mleaf = inf.replace("MLEAF=","")
			elif  "MQRankSum=" in inf:
				mqrs = inf.replace("MQRankSum=","")
			elif  "MQ=" in inf:
				mq = inf.replace("MQ=","")
			elif  "QD=" in inf:
				qd = inf.replace("QD=","")
			elif  "ReadPosRankSum=" in inf:
				rprs = inf.replace("ReadPosRankSum=","")
			elif  "SOR=" in inf:
				sor = inf.replace("SOR=","")





		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for i in (range(len(assb))):

			ass = assb[i]
			nuc = row[4].split(",")[i]




			if row[0] in chrs:



				muts += 1
				sumass += float(ass)
				#print(muts)


				if ("ANN=" in row[7]) and ("LOF=" in row[7]) and ("NMD=" in row[7]):
					print("%s %s %s %s ANN_LOF_NMD" % (row[0],row[1],row[3],nuc) )
					anns = info[-3][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;lof+nmd\n" % (row[0],row[1],row[3],nuc,  row[3],nuc,   a.replace("|",";") ))
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof+nmd\n" % (row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof+nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 


					lofs = info[-2][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))

					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))
					#pass


				elif ("ANN=" in row[7]) and ("LOF=" in row[7]):
					print("%s %s %s %s ANN_LOF" % (row[0],row[1],row[3],nuc) )

					anns = info[-2][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;lof\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, a.replace("|",";") ))
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 

					lofs = info[-1][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))
					#pass

				elif ("LOF=" in row[7]) and ("NMD=" in row[7]):
					print("%s %s %s %s LOF_NMD" % (row[0],row[1],row[3],nuc) )

					lofs = info[-2][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))

					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))
					#pass

					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof_nmd\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof_nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 

				elif ("ANN=" in row[7]) and ("NMD=" in row[7]):
					print("%s %s %s %s ANN_NMD" % (row[0],row[1],row[3],nuc) )
					anns = info[-2][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;nmd\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, a.replace("|",";") ))
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;nmd\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 

					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))
					#pass

				elif "ANN=" in row[7]:
					print("%s %s %s %s ANN" % (row[0],row[1],row[3],nuc) )
					anns = info[-1][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;no\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, a.replace("|",";") ))
					#pass

					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;no\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;no\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 

				elif "LOF=" in row[7]:
					print("%s %s %s %s LOF" % (row[0],row[1],row[3],nuc) )
					lofs = info[-1][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))
					#pass
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 

				elif "NMD=" in row[7]:
					print("%s %s %s %s NMD" % (row[0],row[1],row[3],nuc) )
					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],nuc, a.replace("|",";") ))
				#pass
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_nmd\n" % (row[0],row[1],row[3],nuc,  row[3],nuc, ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],nuc,  row[3],nuc,  ass ,   \
																											an, bqrs, dp, fs, eh, ic, mleac, mleaf, mq, mqrs, qd, rprs, sor,  \
																											mostwanted(anns)[0].replace("|",";"), mostwanted(anns)[1] )) 

	if muts != 0:
		reps5.write("%s;%s;%s\n" % (fi,muts,float(sumass)/muts))
	else:
		reps.write("%s;%s;%s\n" % ( fi,muts,"none" ))    

