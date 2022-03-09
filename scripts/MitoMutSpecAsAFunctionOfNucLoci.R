

rm(list=ls(all=TRUE))

### extra ideas: mtDNA ploidy (fraction of reads covered mt to nuc); deletions (probabilities)
### WE SHOULD PUBLISH EVEN NEGATIVE RESULTS

####################################################################################
###### read mice genotypes from consortium and find THE MUTAH
####################################################################################

Genotypes = read.table('../data/MiceGenotypes/BXD.geno', header=TRUE)

### from Kelley Harris and the paper: C>A fraction was rs52263933 according to our QTL analysis, and haplotypes at this marker should delineate the strains with/without Mutyh mutations as well.

TheMutyhGenotype = Genotypes[Genotypes$Locus == 'rs52263933',]
table(TheMutyhGenotype)
VecOfNamesGenotypes = colnames(TheMutyhGenotype)
TheMutyhGenotypeDataFrame = as.data.frame(t(TheMutyhGenotype))
TheMutyhGenotypeDataFrame$Name = row.names(TheMutyhGenotypeDataFrame)
names(TheMutyhGenotypeDataFrame)=c('NucLoci','NameLoci')

####################################################################################
###### read mtDNA variants and add TheMutyhGenotypeDataFrame column
####################################################################################

MtDna = read.table('../data/all_annmax_single.csv', sep = ';', header=TRUE)
VecOfNamesMtDna = unique(MtDna$Name)

MtDna$SNP = gsub(">",'',MtDna$SNP)

length(intersect(VecOfNamesGenotypes,VecOfNamesMtDna)) # 39 - TOO FEW. ACHTUNG!!!!! need to rename names in MtDna
MtDna$NameNew = MtDna$Name
MtDna$NameNew = gsub("\\_(.*)",'',MtDna$NameNew) #  BXD71_RwwJ > BXD71
MtDna$NameNew = gsub("BXD0","BXD",MtDna$NameNew) #  BXD009 > BXD09
MtDna$NameNew = gsub("BXD0","BXD",MtDna$NameNew) #  BXD09 > BXD9
MtDna$NameNew = gsub("BXD0","BXD",MtDna$NameNew) #  BXD09 > BXD9

VecOfNamesMtDna =sort(unique(MtDna$NameNew))
length(intersect(VecOfNamesGenotypes,VecOfNamesMtDna)) # 94 MUCH BETTER, BUT WE HAVE TO BE SURE IN EACH LINE/NAME 
length(setdiff(VecOfNamesGenotypes,VecOfNamesMtDna)) # 146
# length(setdiff(VecOfNamesMtDna,VecOfNamesGenotypes)) # 6: "BXD221" "BXD222" "BXD224" "BXD227" "C57BL"  "DBA" - ACHTUNG, ancestors

dim(MtDna)
MtDna = merge(MtDna,TheMutyhGenotypeDataFrame, by.x = 'NameNew', by.y = 'NameLoci', all.x = TRUE)
dim(MtDna)
table(MtDna$NucLoci) # who are "H" ACHTUNG!!!! heterozygous???

####################################################################################
###### filtration steps to get normal expected MutSpec
###### each variant can be: 
# 1) difference of ancestal line from RefSeq (should be very common); 
# 2) inherited heteroplasmy (should be more similar in close inbred lines ~ recently diverged lines from the same epoch); 
# 3) de novo (recent) mutation
# 4) noise (sequence artefacts)
####################################################################################

#### filtration of variants (noise, frequency, annotation) ACHTUNG - for now we forget about this.. BUT WE HAVE TO DO IT
MtDna$RefFR = gsub("\\|(.*)",'',MtDna$Allele_specific_forward_reverse_read_counts)
MtDna$AltFR = gsub("(.*)\\|",'',MtDna$Allele_specific_forward_reverse_read_counts)
MtDna$RefF = as.numeric(gsub("\\,(.*)",'',MtDna$RefFR))
MtDna$RefR = as.numeric(gsub("(.*)\\,",'',MtDna$RefFR))
MtDna$AltF = as.numeric(gsub("\\,(.*)",'',MtDna$AltFR))
MtDna$AltR = as.numeric(gsub("(.*)\\,",'',MtDna$AltFR))

## remove low quality (ACHTUNG: epoch):
MtDna = MtDna[(MtDna$AltF + MtDna$AltR) >= 10 & MtDna$AltR > 4 & MtDna$AltF > 4 & (MtDna$RefF + MtDna$RefR) >= 10 & MtDna$RefR > 4 & MtDna$RefF > 4,]

## remove too common on heteroplasmy level (high level of heteroplasmy):
MtDna$MAF = (MtDna$AltF + MtDna$AltR) / (MtDna$AltF + MtDna$AltR + MtDna$RefF + MtDna$RefR)
summary(MtDna$MAF)               # super rare variants! majority are still noise. 
# check methods in nature cancers: https://sciwheel.com/work/item/8185534/resources/8083115/pdf
# what does it mean "Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing" and can we trust it more than VAF?
# MtDna = MtDna[MtDna$MAF < 0.05,] # ACHTUNG: keep high MAF and revert the mutations

summary(MtDna$Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing)
cor.test(MtDna$MAF,MtDna$Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing) # very positive! 
MtDna = MtDna[MtDna$Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing > median(MtDna$Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing),] # !!!

## remove too common within the mouse dataset (changes from RefSeq)

####################################################################################
###### prepare Mutational data: one raw == one inbred line with all mutational data: numbers, frequencies
####################################################################################

#### make a matrix of substitutions
VecOfSubst = c('AT','AG','AC','TA','TG','TC','GT','GA','GC','CT','CG','CA')

FinalMatrix = data.frame(matrix(0,nrow(MtDna),12)); 
names(FinalMatrix) = VecOfSubst

for (column in 1:length(VecOfSubst))
{ # column = 1
subst = VecOfSubst[column]
for (raw in 1:nrow(MtDna))
{ # raw = 1
  if (MtDna$SNP[raw] == subst)  {FinalMatrix[raw,column] = 1}
}
}

FinalMatrix$Total = apply(FinalMatrix[,1:12], 1, FUN = sum)
summary(FinalMatrix$Total) ### good! 

MtDna = cbind(MtDna,FinalMatrix)
names(MtDna)
agg = aggregate(MtDna[46:58], by = list(MtDna$NameNew,MtDna$NucLoci), FUN = sum)
names(agg)[1:2]=c('Name','NucLoci')                
agg$TsTv = (agg$AG + agg$GA + agg$CT + agg$TC)/ (agg$AT + agg$TA + agg$AC + agg$CA + agg$GC + agg$CG + agg$GT + agg$TG)
agg$FrAG = agg$AG/agg$Total
agg$FrGA = agg$GA/agg$Total
agg$FrCT = agg$CT/agg$Total
agg$FrTC = agg$TC/agg$Total
agg$FrGT = agg$GT/agg$Total
agg$FrCA = agg$CA/agg$Total

summary(agg)

####################################################################################
###### GWAS EXPECTED: an association between NucLoci and MutSpec (number and spectrum)
####################################################################################

ForGwas = agg[agg$NucLoci %in% c('B','D'),] # in D we expect higher oxidative damage
names(ForGwas)
agg1 = aggregate(ForGwas[3:22], by = list(ForGwas$NucLoci), FUN = mean)  # ACHTUNG: doesn't look normal mtDNA MutSpec
# we don't go down before we have nice MutSpec

wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$Total,ForGwas[ForGwas$NucLoci == 'D',]$Total)
wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$FrAG,ForGwas[ForGwas$NucLoci == 'D',]$FrAG)
wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$FrTC,ForGwas[ForGwas$NucLoci == 'D',]$FrTC)
wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$FrCT,ForGwas[ForGwas$NucLoci == 'D',]$FrCT)
wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$FrGA,ForGwas[ForGwas$NucLoci == 'D',]$FrGA)

wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$FrGT,ForGwas[ForGwas$NucLoci == 'D',]$FrGT)
wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$FrCA,ForGwas[ForGwas$NucLoci == 'D',]$FrCA)

t.test(ForGwas[ForGwas$NucLoci == 'B',]$FrGT + ForGwas[ForGwas$NucLoci == 'B',]$FrCA, ForGwas[ForGwas$NucLoci == 'D',]$FrGT + ForGwas[ForGwas$NucLoci == 'D',]$FrCA)
wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$FrGT + ForGwas[ForGwas$NucLoci == 'B',]$FrCA, ForGwas[ForGwas$NucLoci == 'D',]$FrGT + ForGwas[ForGwas$NucLoci == 'D',]$FrCA)

wilcox.test(ForGwas[ForGwas$NucLoci == 'B',]$TsTv,ForGwas[ForGwas$NucLoci == 'D',]$TsTv)

# A bit higher GT an CA in D vs B
summary(ForGwas[ForGwas$NucLoci == 'D',]$FrGT)
summary(ForGwas[ForGwas$NucLoci == 'B',]$FrGT)

summary(ForGwas[ForGwas$NucLoci == 'D',]$FrCA)
summary(ForGwas[ForGwas$NucLoci == 'B',]$FrCA)

# A bit higher Total in D vs B
summary(ForGwas[ForGwas$NucLoci == 'D',]$Total)
summary(ForGwas[ForGwas$NucLoci == 'B',]$Total)

# also a bit higher TC (Ah>Gh) in D vs B
summary(ForGwas[ForGwas$NucLoci == 'D',]$FrTC)
summary(ForGwas[ForGwas$NucLoci == 'B',]$FrTC)

# also a bit higher TC (Ah>Gh) in D vs B
summary(ForGwas[ForGwas$NucLoci == 'D',]$FrGA)
summary(ForGwas[ForGwas$NucLoci == 'B',]$FrGA)
