import time
import timeit
import os







pgen_path = '_DATA FROM KELLY&DAVID\\BXD_Sequencing_Summary_Statistics_21_Jan_2019.csv'
pgen = open(pgen_path,'r')




path3="new.csv"
reps3 = open(path3, "w")


pgen.readline().strip()
reps3.write("Index_Row;HudsonAlpha ID;Alternate ID;Strain;Generation_JAX;DOB;Age at death (days);Parents ;Type;Source;Team;Organism ;Sex ;Tissue;longranger_version;instrument_id; gems_detected ; mean_dna_per_gem ;bc_on_whitelist;bc_mean_qscore;n50_linked_reads_per_molecule;corrected_loaded_mass_ng; snps_phased ;genes_phased_lt_100kb; longest_phase_block ; n50_phase_block ; molecule_length_mean ; molecule_length_stddev ; number_reads ;median_insert_size; mean_depth ;zero_coverage;mapped_reads;pcr_duplication;on_target_bases;r1_q20_bases_fract;r1_q30_bases_fract;r2_q20_bases_fract;r2_q30_bases_fract;si_q20_bases_fract;si_q30_bases_fract;bc_q20_bases_fract;bc_q30_bases_fract;large_sv_calls;short_deletion_calls;Expanded name;Official name;GeneNetwork name;Year breeding started;Date production began at JAX;JAX Stock No (RRID);Availability (as of March 2019);Availability (as of October 2019);Epoch;Method of derivation;Mitocondrial origin (if known);October 2018 reported generation at JAX ;March 2019 reported generation at JAX;October 2019 reported generation at JAX;Febuary 2019 generation at UTHSC ;October 2019 generation at UTHSC;Genotyped on Illumina 13377 SNP array;Genotyped on Affymetrix MouseDiversityArray600K;Generation at Affymetrix MouseDiversityArray600K;Generation at MUGA array genotyping;Generation at MegaMUGA array genotyping;Generation at GigaMUGA array genotyping;Generation went extinct (if known);Backcross and RIX-derived breeding notes;Any genotypes available\n")


while True:

    read = pgen.readline().strip()
    if read == "":
        break
        
    name = read.split(';')[3]

    filepath = '_DATA FROM KELLY&DAVID\\Table_S1_20_Nov_2020.csv'
    file = open(filepath,'r')


    mark = 0

    while True:
        raw = file.readline().strip()
        if raw == '':
            break

        name2 = raw.split(';')[1]

        if name2.replace('/TyJ','').replace('/Ty','').replace('/RwwJ','').replace('/2RwwJ','') == name or "0"+name2.replace('BXD','').replace('/TyJ','').replace('/Ty','').replace('/RwwJ','').replace('/2RwwJ','') == name.replace('BXD','')or "00"+name2.replace('BXD','').replace('/TyJ','').replace('/Ty','').replace('/RwwJ','').replace('/2RwwJ','') == name.replace('BXD','')                  :

            reps3.write("%s;%s\n" % (read,raw))

            mark = 1

    if mark == 0:

        reps3.write("%s\n" % (read))