import os 

#files = ["Control-Polg-1","Control-Polg-2","Control-Polg-3","Parkin-Polg-1","Parkin-Polg-2","Parkin-Polg-3"]
#files = ["MGI-1","MGI-2","MGI-3","MGI-4","MGI-97","MGI-98","MGI-99","MGI-100","MGI-101","MGI-102","MGI-103","MGI-104"]
files=os.listdir("data\\M\\")

#prepath = "C:\\Users\\Abathur\\Desktop\\ParkPolGbam\\mutect_first\\"
prepath = "C:\\Users\\Abathur\\Desktop\\Science\\InbredMice\\data\\M\\"
prepath2 = "C:\\Users\\Abathur\\Desktop\\Science\\InbredMice\\data\\M_split\\"

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
reps6.write("File;Name;Chr;Pos;Ref;Alt;SNP;Allele_specific_forward_reverse_read_counts;variant_allele_frequency;Approximate_read_depth;Number_of_events_in_this_haplotype;median_base_quality;median_fragment_length;median_mapping_quality;median_distance_from_end_of_read;Number_of_alt_reads_original_alignment_doesnt_match_the_current_contig;negative_log10_population_allele_frequencies_of_alt_alleles;Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing;Allele;Annotation;Annotation_Impact;Gene_Name;Gene_ID;Feature_Type;Feature_ID;Transcript_BioType;Rank;HGVS.c;HGVS.p;cDNA.pos/cDNA.length;CDS.pos/CDS.length;AA.pos/AA.length;Distance;ERRORS/WARNINGS/INFO;Score;lof\n") 
	




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
	reps4.write("Chr;Pos;Ref;Alt;SNP;Allele_specific_forward_reverse_read_counts;variant_allele_frequency;Approximate_read_depth;Number_of_events_in_this_haplotype;median_base_quality;median_fragment_length;median_mapping_quality;median_distance_from_end_of_read;Number_of_alt_reads_original_alignment_doesnt_match_the_current_contig;negative_log10_population_allele_frequencies_of_alt_alleles;Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing;Allele;Annotation;Annotation_Impact;Gene_Name;Gene_ID;Feature_Type;Feature_ID;Transcript_BioType;Rank;HGVS.c;HGVS.p;cDNA.pos/cDNA.length;CDS.pos/CDS.length;AA.pos/AA.length;Distance;ERRORS/WARNINGS/INFO;Score;lof\n") 
	

	
	name = "_".join(fi.split("_phased_")[0].split("_")[1:])



	##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.">
	##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description="Number of reads with unique start and mate end positions for each alt at a variant site">
	##INFO=<ID=CONTQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to contamination">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
	##INFO=<ID=ECNT,Number=1,Type=Integer,Description="Number of events in this haplotype">
	##INFO=<ID=GERMQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not germline variants">
	##INFO=<ID=MBQ,Number=R,Type=Integer,Description="median base quality">
	##INFO=<ID=MFRL,Number=R,Type=Integer,Description="median fragment length">
	##INFO=<ID=MMQ,Number=R,Type=Integer,Description="median mapping quality">
	##INFO=<ID=MPOS,Number=A,Type=Integer,Description="median distance from end of read">
	##INFO=<ID=NALOD,Number=A,Type=Float,Description="Negative log 10 odds of artifact in normal with same allele fraction as tumor">
	##INFO=<ID=NCount,Number=1,Type=Integer,Description="Count of N bases in the pileup">
	##INFO=<ID=NLOD,Number=A,Type=Float,Description="Normal log 10 likelihood ratio of diploid het or hom alt genotypes">
	##INFO=<ID=OCM,Number=1,Type=Integer,Description="Number of alt reads whose original alignment doesn't match the current contig.">
	##INFO=<ID=PON,Number=0,Type=Flag,Description="site found in panel of normals">
	##INFO=<ID=POPAF,Number=A,Type=Float,Description="negative log 10 population allele frequencies of alt alleles">
	##INFO=<ID=ROQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to read orientation artifact">
	##INFO=<ID=RPA,Number=R,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
	##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
	##INFO=<ID=SEQQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not sequencing errors">
	##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
	##INFO=<ID=STRANDQ,Number=1,Type=Integer,Description="Phred-scaled quality of strand bias artifact">
	##INFO=<ID=STRQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors">
	##INFO=<ID=TLOD,Number=A,Type=Float,Description="Log 10 likelihood ratio score of variant existing versus not existing">






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
		assb = info[0].replace("AS_SB_TABLE=","").split("|") # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   1665,1643|9,4|3,0
		allas = 0
		vaf = [] # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		# print(fi)
		# print(info) 

		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    считаем все риды всех аллелей (forward и reverse) 
		for ass in assb:
			ab = ass.split(",")
			sumab = int(ab[0]) + int(ab[1])
			allas += sumab


		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    считаем отношение ридов данной аллели (forward и reverse) ко всем ридам 
		for ass in assb:
			ab = ass.split(",")
			sumab = int(ab[0]) + int(ab[1])
			if allas != 0:
				vaf.append(str(float(float(sumab)/float(allas))))
			else:
				vaf.append(str(0))
		


		for i in (range(len(assb[1:]))): # для всех альтернативных аллелей берем их нуклеотид и посчитанную частоту
			ass = assb[i+1]
			vafass = vaf[i+1]  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			nuc = row[4].split(",")[i]




			if row[0] in chrs:



				muts += 1
				sumass += float(vafass)
				#print(muts)


				if ("ANN=" in row[7]) and ("LOF=" in row[7]) and ("NMD=" in row[7]):
					print("%s %s %s %s ANN_LOF_NMD" % (row[0],row[1],row[3],row[4]) )
					anns = info[-3][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;lof+nmd\n" % (row[0],row[1],row[3],row[4],  row[3],row[4],   a.replace("|",";") ))
					# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof+nmd\n" % (row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof+nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 

					lofs = info[-2][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))

					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))
					#pass


				elif ("ANN=" in row[7]) and ("LOF=" in row[7]):
					print("%s %s %s %s ANN_LOF" % (row[0],row[1],row[3],row[4]) )

					anns = info[-2][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;lof\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], a.replace("|",";") ))
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass)  ,info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""),info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""),info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""),  mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;lof\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 

					lofs = info[-1][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))
					#pass

				elif ("LOF=" in row[7]) and ("NMD=" in row[7]):
					print("%s %s %s %s LOF_NMD" % (row[0],row[1],row[3],row[4]) )

					lofs = info[-2][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))

					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))
					#pass

					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof_nmd\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass)  ,info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""),info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""),info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""),  mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof_nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 

				elif ("ANN=" in row[7]) and ("NMD=" in row[7]):
					print("%s %s %s %s ANN_NMD" % (row[0],row[1],row[3],row[4]) )
					anns = info[-2][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;nmd\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], a.replace("|",";") ))
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;nmd\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass)  ,info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""),info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""),info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""),  mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 

					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))
					#pass

				elif "ANN=" in row[7]:
					print("%s %s %s %s ANN" % (row[0],row[1],row[3],row[4]) )
					anns = info[-1][4:].split(",")
					for a in anns:
						#print("    %s" % a.split("|"))
						reps.write("%s;%s;%s;%s;%s>%s;%s;no\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], a.replace("|",";") ))
					#pass
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;no\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass)  ,info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""),info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""),info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""),  mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;no\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 

				elif "LOF=" in row[7]:
					print("%s %s %s %s LOF" % (row[0],row[1],row[3],row[4]) )
					lofs = info[-1][5:-1].split(",")
					for a in lofs:
						#print("    %s" % a.split("|"))
						reps2.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))
					#pass
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass)  ,info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""),info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""),info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""),  mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_lof\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 

				elif "NMD=" in row[7]:
					print("%s %s %s %s NMD" % (row[0],row[1],row[3],row[4]) )
					nmds = info[-1][5:-1].split(",")
					for a in nmds:
						#print("    %s" % a.split("|"))
						reps3.write("%s;%s;%s;%s;%s\n" % (row[0],row[1],row[3],row[4], a.replace("|",";") ))
				#pass
					reps4.write("%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_nmd\n" % (row[0],row[1],row[3],row[4],  row[3],row[4], '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass)  ,info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""),info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""),info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""),  mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] ))
					reps6.write("%s;%s;%s;%s;%s;%s;%s>%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;only_nmd\n" % (fi.split(".pacsortmark")[0], name, row[0],row[1],row[3],row[4],  row[3],row[4],  '%s|%s' % (assb[0],ass) ,   '%s|%s' % (vaf[0],vafass) , \
																										info[1].replace("DP=",""),info[2].replace("ECNT=",""),info[3].replace("MBQ=",""), \
																										info[4].replace("MFRL=",""),info[5].replace("MMQ=",""),info[6].replace("MPOS=",""), \
																										info[7].replace("OCM=",""),info[8].replace("POPAF=",""),info[9].replace("TLOD=",""), \
																										mostwanted(anns)[0].replace("|",";"),mostwanted(anns)[1] )) 
	
	if muts != 0:
		reps5.write("%s;%s;%s\n" % (fi,muts,float(sumass)/muts))
	else:
		reps.write("%s;%s;%s\n" % ( fi,muts,"none" ))    

