import time
import timeit
import os







pgen_path = 'all_annmax_single.csv'
pgen = open(pgen_path,'r')




path3="merged.csv"
reps3 = open(path3, "w")


pgen.readline().strip()
reps3.write("File;Name;Chr;Pos;Ref;Alt;SNP;Allele_specific_forward_reverse_read_counts;variant_allele_frequency;Approximate_read_depth;Number_of_events_in_this_haplotype;median_base_quality;median_fragment_length;median_mapping_quality;median_distance_from_end_of_read;Number_of_alt_reads_original_alignment_doesnt_match_the_current_contig;negative_log10_population_allele_frequencies_of_alt_alleles;Log_10_likelihood_ratio_score_of_variant_existing_versus_not_existing;Allele;Annotation;Annotation_Impact;Gene_Name;Gene_ID;Feature_Type;Feature_ID;Transcript_BioType;Rank;HGVS.c;HGVS.p;cDNA.pos/cDNA.length;CDS.pos/CDS.length;AA.pos/AA.length;Distance;ERRORS/WARNINGS/INFO;Score;lof;OrdId;IdNum;File;Name;IdNum2;Index_Row;HudsonAlpha ID;Alternate ID;Strain;Generation_JAX;DOB;Age at death (days);Parents ;Type;Source;Team;Organism ;Sex ;Tissue;longranger_version;instrument_id; gems_detected ; mean_dna_per_gem ;bc_on_whitelist;bc_mean_qscore;n50_linked_reads_per_molecule;corrected_loaded_mass_ng; snps_phased ;genes_phased_lt_100kb; longest_phase_block ; n50_phase_block ; molecule_length_mean ; molecule_length_stddev ; number_reads ;median_insert_size; mean_depth ;zero_coverage;mapped_reads;pcr_duplication;on_target_bases;r1_q20_bases_fract;r1_q30_bases_fract;r2_q20_bases_fract;r2_q30_bases_fract;si_q20_bases_fract;si_q30_bases_fract;bc_q20_bases_fract;bc_q30_bases_fract;large_sv_calls;short_deletion_calls;Expanded name;Official name;GeneNetwork name;Year breeding started;Date production began at JAX;JAX Stock No (RRID);Availability (as of March 2019);Availability (as of October 2019);Epoch;Method of derivation;Mitocondrial origin (if known);October 2018 reported generation at JAX ;March 2019 reported generation at JAX;October 2019 reported generation at JAX;Febuary 2019 generation at UTHSC ;October 2019 generation at UTHSC;Genotyped on Illumina 13377 SNP array;Genotyped on Affymetrix MouseDiversityArray600K;Generation at Affymetrix MouseDiversityArray600K;Generation at MUGA array genotyping;Generation at MegaMUGA array genotyping;Generation at GigaMUGA array genotyping;Generation went extinct (if known);Backcross and RIX-derived breeding notes;Any genotypes available\n")


while True:

    read = pgen.readline().strip()
    if read == "":
        break
        
    name = read.split(';')[0]

    filepath = 'new.csv'
    file = open(filepath,'r')


    mark = 0

    while True:
        raw = file.readline().strip()
        if raw == '':
            break

        name2 = raw.split(';')[2]

        if name2.strip() == name.strip():

            reps3.write("%s;%s\n" % (read,raw))

            mark = 1

    if mark == 0:

        reps3.write("%s\n" % (read))