import os 
import shutil


prepath = "/mnt/c/Users/Abathur/Desktop/Science/InbredMice/new/Mutect/"
   


files=os.listdir(prepath)

refseq = "/mnt/c/Users/Abathur/Desktop/Science/_RefSeqs/GRCm38_ChrM/Mus_musculus.GRCm38.dna.primary_assembly_ChrM.fa"



for fi in files:


	fsrc = "%s" % (prepath)
	fdst = "/mnt/c/Users/Abathur/Desktop/Science/InbredMice/new/temp3/"



	os.system('tools/vt/./vt  decompose -s %s%s -o %s%s' % (fsrc, fi, fdst, fi.replace(".vcf",".decomposed.vcf"))) 
	os.system('tools/vt/./vt  normalize %s%s -r %s | tools/vt/./vt uniq - -o %s%s' % (fdst, fi.replace(".vcf",".decomposed.vcf"), refseq, fdst, fi.replace(".vcf",".decomposed.normalized.uniq.vcf")))