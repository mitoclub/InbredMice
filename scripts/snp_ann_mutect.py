import os 
import shutil


#prepath = "C:\\Users\\Abathur\\Desktop\\ParkPolGbam\\first_snp\\"
#prepath = "C:\\Users\\Abathur\\Desktop\\ParkPolGbam\\second_snp\\"
# prepath = "/home/jester/data/temp3/"

prepath = "/home/jester/_Science/new/temp3/"








files=os.listdir("temp3")





for fi in files:


	# GENERATE FILE avida.cfg

	fsrc = "%s%s" % (prepath, fi)
	fdst = "/home/jester/_Science/new/temp4/%s" % (fi)
	# shutil.copyfile(fsrc, fdst)

	with open (fsrc, 'r') as f:
		old_data = f.read()

	new_data = old_data.replace('chrM', 'MT')
	# new_data = new_data.replace('WORLD_Y 60', 'WORLD_Y %s' % m)

	with open (fdst, 'w') as f:
		f.write(new_data)





	os.system('java -Xmx16g -jar /home/jester/tools/snpEff/snpEff.jar -no-upstream -no-downstream -v GRCm38.99 %s  > /home/jester/_Science/new/MutectVTSnpEff/%s' % (fdst, fi))
	