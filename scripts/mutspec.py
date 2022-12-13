import time
import timeit
import os





filepath = 'new.csv'
file = open(filepath,'r')




subs = ["A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>T", "C>G"]

# T>A    1160
# C>T    1072
# T>C    1072
# A>T    1040
# C>A    1026
# A>C     959
# A>G     904
# C>G     728
# G>A     657
# T>G     634
# G>T     589
# G>C     201





path3="mutspec.csv"
reps3 = open(path3, "w")


file.readline().strip()
reps3.write("OrdId;IdNum;File;Name;IdNum2;Index_Row;A>T;A>G;A>C;T>A;T>G;T>C;G>A;G>T;G>C;C>A;C>T;C>G;HudsonAlpha ID;Alternate ID;Strain;Generation_JAX;DOB;Age at death (days);Parents ;Type;Source;Team;Organism ;Sex ;Tissue;longranger_version;instrument_id; gems_detected ; mean_dna_per_gem ;bc_on_whitelist;bc_mean_qscore;n50_linked_reads_per_molecule;corrected_loaded_mass_ng; snps_phased ;genes_phased_lt_100kb; longest_phase_block ; n50_phase_block ; molecule_length_mean ; molecule_length_stddev ; number_reads ;median_insert_size; mean_depth ;zero_coverage;mapped_reads;pcr_duplication;on_target_bases;r1_q20_bases_fract;r1_q30_bases_fract;r2_q20_bases_fract;r2_q30_bases_fract;si_q20_bases_fract;si_q30_bases_fract;bc_q20_bases_fract;bc_q30_bases_fract;large_sv_calls;short_deletion_calls;Expanded name;Official name;GeneNetwork name;Year breeding started;Date production began at JAX;JAX Stock No (RRID);Availability (as of March 2019);Availability (as of October 2019);Epoch;Method of derivation;Mitocondrial origin (if known);October 2018 reported generation at JAX ;March 2019 reported generation at JAX;October 2019 reported generation at JAX;Febuary 2019 generation at UTHSC ;October 2019 generation at UTHSC;Genotyped on Illumina 13377 SNP array;Genotyped on Affymetrix MouseDiversityArray600K;Generation at Affymetrix MouseDiversityArray600K;Generation at MUGA array genotyping;Generation at MegaMUGA array genotyping;Generation at GigaMUGA array genotyping;Generation went extinct (if known);Backcross and RIX-derived breeding notes;Any genotypes available\n")


while True:

    read = file.readline().strip()
    if read == "":
        break
        
    name = read.split(';')[2]



    

    reps3.write("%s;" % (';'.join(read.split(';')[:6])))

    for sub in subs:


        pgen_path = 'all_annmax_single.csv'
        pgen = open(pgen_path,'r')

        count = 0

        while True:
            raw = pgen.readline().strip()
            if raw == '':
                break

            name2 = raw.split(';')[0]

            if name2.strip() == name.strip():



                if raw.split(';')[6] == sub:

                    count+=1




        reps3.write("%s;" % (count))


    reps3.write("%s" % (';'.join(read.split(';')[6:])))

    reps3.write("\n")       