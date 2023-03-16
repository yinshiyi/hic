####
# use this to mine existing data
aws s3 ls s3://egenesis-data-raw --recursive | grep HiC | cut -c 32- | head


##############
# using this to check if the reads are paired
fastq_pair <(zcat sus_scrofa_ovary_ref_bestcontigs_S3HiC_R1.fastq.gz) <(zcat sus_scrofa_ovary_ref_bestcontigs_S3HiC_R2.fastq.gz)

##############
cut -f1,2,3,4,5,6,7,8,9,10 ~/mnt/pigliveroutput/hic_results/data/sample1/sample1.allValidPairs | awk '{OFS = "\t";print $1,$4,$2,$3,$9,$7,$5,$6,$10,$8}' | awk '{ OFS = "\t"; if($2=="+") $2=0 ; if($2=="-") $2=16 ;if($6=="+") $6=0 ; if($6=="-") $6=16 ; $1=""; $5=0; $9=0;print $0}' |  awk '{gsub("chr","",$2);print}' |awk '{gsub("chr","",$6);print}' | sort -k2,2d -n -k6,6d -n | cut -f2,3,4,5,6,7,8,9,10 > ~/mnt/pigliveroutput/hic_results/data/sample1/reformat.reorder.correct.sorted.pairs 
#the properly sorted file did not change the missing cis contact map
java -Xmx2g -jar ~/myJuicerDir/scripts/common/juicer_tools.jar pre ~/mnt/pigliveroutput/hic_results/data/sample1/reformat.reorder.correct.sorted.pairs pigliver.hic ~/susScr11.chrom.sizes -d true
# the d option will only look at cis contact, that is when the error code happens

##############
# need to build ensembl genome reference
bowtie2-build --threads 32 susScr11.fa ./ss11/Sscrofa11.1
/home/shiyi/mnt/2021genome/bowtie

###########
cat SRR10764659_1.fastq.gz SRR10764658_1.fastq.gz > ./muscle/sample1/muscle_R1.fastq.gz
cat SRR10764659_2.fastq.gz SRR10764658_2.fastq.gz > ./muscle/sample1/muscle_R2.fastq.gz
./HiC-Pro-master/bin/HiC-Pro -i ~/mnt/muscle/ -o ~/mnt/pigmuscleoutput/ -c ~/mnt/configMboI.txt 

###############
./HiC-Pro-master/bin/utils/digest_genome.py -r A^AGCTT -o pig_hindiii.bed <(zcat ~/susScr11.fa
.gz)
./HiC-Pro-master/bin/utils/digest_genome.py -r dpnii -o pig_DpnII.bed susScr11.fa

############33
# downsampling
seqtk sample -s 123 liver/sample1/ERR2652854_1.fastq.gz 100000 |gzip > ./liver2/sample2/ERR265284_short1.fastq.gz
seqtk sample -s100 read1.fq 10000 > sub1.fq
############
# merging the two sequencing runs analysis
# Based on what I gathered here and what the Shuhong Zhao wrote https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4256199
# After mapping, the singleton, multi mapped, dumped, dangling and self-circle paired-end reads and PCR duplication were all removed by HiC-Pro with the parameter “-s proc_hic”. 
# Next, files with.bwt2pairs.validPairs suffix from different libraries were moved to a new folder as the input data for subsequent analysis. 
# And then merging multiple libraries, building raw inter-/intra-chromosomal contact maps and running ICE normalization on contact maps were performed using the programs in order of HiC-Pro parameter “-s merge_persample”, “-s build_contact_maps” and “-s ice_norm”.
# regarding merging samples
# https://groups.google.com/g/hic-pro/c/tRVfjs8eigs
cat ~/mnt/pigmuscleoutput_secondhalf/hic_results/data/sample1/sample1.allValidPairs \
    ~/mnt/pigmuscleoutput_firsthalf/hic_results/data/sample1/sample1.allValidPairs > pigmuscle.allvalidpairs
cut -f1,2,3,4,5,6,7,8,9,10 ~/mnt/pigmuscle.allvalidpairs | awk '{OFS = "\t";print $1,$4,$2,$3,$9,$7,$5,$6,$10,$8}' | awk '{ OFS = "\t"; if($2=="+") $2=0 ; if($2=="-") $2=16 ;if($6=="+") $6=0 ; if($6=="-") $6=16 ; $1=""; $5=0; $9=1;print $0}' | awk '{gsub("chr","",$2);print}' |awk '{gsub("chr","",$6);print}' | sort -k2,2d -n -k6,6d -n |  awk '{OFS = "\t";print $1,$2,$3,$4,$5,$6,$7,$8}' > ~/mnt/pigmuscle.reformat.reorder.correct.sorted.pairs

# juicer
java -Xmx10g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar pre ~/mnt/pigmuscle.reformat.reorder.correct.sorted.pairs ~/mnt/pig_muscle_full.hic ~/mnt/susScr11.chrom.sizes


#######
# hicuup
# https://github.com/aidenlab/juicer/wiki/HiCCUPS
java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups
java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups \
--cpu --threads 4 \
-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 -c 6 \
--ignore_sparsity /home/shiyi/mnt/pig_muscle_full.hic \
pig.hic.loops
##################
java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups \
--cpu --threads 4 \
-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 -c 17 \
--ignore_sparsity /home/shiyi/mnt/pig_muscle_full.hic \
GSH18_pig.hic.loops

java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups \
--cpu --threads 4 \
-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 \
--ignore_sparsity /home/shiyi/mnt/pig_muscle_full.hic \
full_pig.hic.loops


genomeID="Sscrofa11.1"
hic_file_path="/home/shiyi/mnt/pig_muscle_full.hic"
juicer_tools_path="/home/shiyi/mnt/myJuicerDir/scripts/common/juicer_tools"
bed_file_dir="/home/shiyi/mnt/eGen_HiC/bedfile/pig_MboI.bed"

############
cat ./CTCFmining/laddingpad.gtf | gtf2bed | awk '{ OFS = "\t"; $5=0 ;print $0}' | head

# generate bed file for visulization
awk '{ OFS = "\t"; if($1=="chr6"&&$4>=58793153&&$5<=58902052&&$3=="transcript") print $0}' susScr11.ncbiRefSeq.gtf | gtf2bed | awk '{ OFS = "\t"; $5=0 ;print $0}' | cut -f1,2,3,4,5,6 > landingpadgenes.bed
Rscript bedfile.R landingpadgenes.bed landingpadgenes_filtered.bed
# 6:111,433,875-111,562,238
#    awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' |
# a solution from the internet
# histroical debug
head -n 118720 /home/shiyi/mnt/2021genome/2021payload_inserted_ss11/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.gtf |    tail +118720 | awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' | gtf2bed

# these code using gtf2bed does not work well with ensembl
cat /home/shiyi/mnt/2021genome/2021payload_inserted_ss11/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.gtf | \
    awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' | \
    gtf2bed | awk '{ OFS = "\t"; $5=0 ;print $0}' | \
    awk '{ OFS = "\t"; if($1=="6"&&$2>=111433875&&$3<=111562238) print $0}' | \
    cut -f1,2,3,4,5,6 > legacy_genome_landingpadgenes.bed
# using bob's rscript to deal with ensembl, to do list
cd /home/shiyi/mnt/gtf2bed
perl gtf2bed2.pl /home/shiyi/mnt/2021genome/2021payload_inserted_ss11/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.gtf genome_gene.bed
# a modified perl script will make the gene name appear
awk '{ OFS = "\t"; if($1=="6"&&$2>=111003875&&$3<=111582238) print $0}' genome_gene.bed > legacy_genome_landingpadgenes.bed
/home/shiyi/mnt/2021genome/ss11_illustor.bed
awk '{ OFS = "\t"; if($1=="12"&&$2>=40789561&&$3<=40808560) print $0}' genome_gene.bed > ccl2_landingpadgenes.bed
awk '{ OFS = "\t"; if($1=="12"&&$2>=40000061&&$3<=41008560) print $0}' /home/shiyi/mnt/2021genome/ss11_illustor.bed > /home/shiyi/mnt/eGen_HiC/ccl2_landingpadgenes.bed

# 
perl gtf2bed2.pl /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.104.gtf genome_gene_regular.bed
awk '{ OFS = "\t"; if($1=="6"&&$2>=58673465&&$3<=59509568) print $0}' genome_gene_regular.bed > regular_genome_landingpadgenes.bed
# ccl2 
# chr12:40789561-40808560
Rscript bedfile.R landingpadgenes.bed landingpadgenes_filtered.bed
#########
# adding chr6 to fit ucsc/ncbi genome browser
awk 'OFS="\t" {$1="chr"$1; print}' /home/shiyi/mnt/2021genome/illustor.bed > /home/shiyi/mnt/eGen_HiC/landingpad.bed

#######
# hichip data

nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/HiChIP/HiChIP_fastqs/*_R{1,2}_001.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --dnase --min_cis 1000 \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000

# aligning to pl15s genome
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-164_R{1,2}.fastq.gz' \
    --bwa_index /home/shiyi/mnt/2021genome/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.fa    \
    --fasta /home/shiyi/mnt/2021genome/2021payload_inserted_ss11/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.fa \
    --dnase --min_cis 1000 \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference

# output looks great

##################
# merge the fastq file before running hic pro
~/mnt/hic/results/hicpro/valid_pairs
# cat ~/mnt/hic/results/hicpro/valid_pairs/DTG-HiChIP-164_S6_L00?.allValidPairs 
cat /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-164_S6_L00?_R1_001.fastq.gz > DTG-HiChIP-164_R1.fastq.gz
cat /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-164_S6_L00?_R2_001.fastq.gz > DTG-HiChIP-164_R2.fastq.gz
# zcat *R1*.fastq.gz | wc -l

################################
# reformat the dataframe
cut -f1,2,3,4,5,6,7,8,9,10 /home/shiyi/mnt/hic/results/hicpro/valid_pairs/DTG-HiChIP-164_R.allValidPairs | \
awk '{OFS = "\t";print $1,$4,$2,$3,$9,$7,$5,$6,$10,$8}' | \
awk '{ OFS = "\t"; if($2=="+") $2=0 ; if($2=="-") $2=16 ;if($6=="+") $6=0 ; if($6=="-") $6=16 ; $1=""; $5=0; $9=1;print $0}' |  \
sort -k2,2d -k3,3n  | \
awk '{OFS = "\t";print $1,$2,$3,$4,$5,$6,$7,$8}'  > \
/home/shiyi/mnt/hic/results/hicpro/valid_pairs/DTG-HiChIP-164_reformat.reorder.correct.sorted.pairs


#################
# making the hic file
# https://github.com/aidenlab/juicer/wiki/Pre
java -Xmx2g -jar /home/shiyi/mnt/myJuicerDir/scripts/common/juicer_tools.jar \
pre /home/shiyi/mnt/hic/results/hicpro/valid_pairs/DTG-HiChIP-164_reformat.reorder.correct.sorted.pairs \
/home/shiyi/mnt/hic/results/hicpro/DTG-164.hic \
/home/shiyi/mnt/2021genome/2021payload_inserted_ss11/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.chromSize \
-d true
# only calculate intra chromosomes

##################
# plotting the looping

java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups \
--cpu --threads 4 \
-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000  \
--ignore_sparsity /home/shiyi/mnt/hic/results/hicpro/DTG-164.hic \
/home/shiyi/mnt/hic/results/hicpro/DTG-164.loops

# default parameter
# 8 loops consider decrease the threshold
java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups \
--cpu --threads 4 \
--ignore_sparsity /home/shiyi/mnt/hic/results/hicpro/DTG-164.hic \
/home/shiyi/mnt/hic/results/hicpro/DTG-164.loops

###############
cat /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-162_S10_L00?_R1_001.fastq.gz /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-163_S9_L00?_R1_001.fastq.gz > /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-ctcf_R1.fastq.gz
cat /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-162_S10_L00?_R2_001.fastq.gz /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-163_S9_L00?_R2_001.fastq.gz > /home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-ctcf_R2.fastq.gz

#############
cd /home/shiyi/mnt/hic
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/HiChIP/HiChIP_fastqs/DTG-HiChIP-ctcf_R{1,2}.fastq.gz' \
    --bwa_index /home/shiyi/mnt/2021genome/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.fa    \
    --fasta /home/shiyi/mnt/2021genome/2021payload_inserted_ss11/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.fa \
    --dnase --min_cis 1000 \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference

cut -f1,2,3,4,5,6,7,8,9,10 /home/shiyi/mnt/hic/results/hicpro/valid_pairs/DTG-HiChIP-ctcf_R.allValidPairs | \
awk '{OFS = "\t";print $1,$4,$2,$3,$9,$7,$5,$6,$10,$8}' | \
awk '{ OFS = "\t"; if($2=="+") $2=0 ; if($2=="-") $2=16 ;if($6=="+") $6=0 ; if($6=="-") $6=16 ; $1=""; $5=0; $9=1;print $0}' |  \
sort -k2,2d -k3,3n  | \
awk '{OFS = "\t";print $1,$2,$3,$4,$5,$6,$7,$8}'  > \
/home/shiyi/mnt/eGen_HiC/DTG-HiChIP-ctcf_reformat.reorder.correct.sorted.pairs

java -Xmx2g -jar /home/shiyi/mnt/myJuicerDir/scripts/common/juicer_tools.jar \
pre \
/home/shiyi/mnt/eGen_HiC/DTG-HiChIP-ctcf_reformat.reorder.correct.sorted.pairs \
/home/shiyi/mnt/eGen_HiC/DTG-ctcf.hic \
/home/shiyi/mnt/2021genome/2021payload_inserted_ss11/Sus_scrofa11.1_PL15s_integrated_chr6_reversed_other_nonreversed.chromSize \
-d true

java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups \
--cpu --threads 14 \
-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000  \
--ignore_sparsity \
/home/shiyi/mnt/eGen_HiC/DTG-ctcf.hic \
/home/shiyi/mnt/eGen_HiC/DTG-ctcf.loops



###############
#WT
cd /home/shiyi/mnt/hic
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_R{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --restriction_site 'G^ANTC ^GATC C^TNAG T^TAA' \
    --ligation_site 'ANTCGANT,ANTCGATC,ANTCCTNA,ANTCTTA,GATCGANT,GATCGATC,GATCCTNA,GATCTTA,TNAGGANT,TNAGGATC,TNAGCTNA,TNAGTTA,TAAGANT,TAAGATC,TAACTNA,TAATTA' \
    --restriction_fragments /home/shiyi/mnt/eGen_HiC/phasegenomics/pig_pg.bed \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney \
    -resume curious_pare

############3
# use the new bed file
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_R{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --restriction_site 'G^ANTC ^GATC C^TNAG T^TAA' \
    --ligation_site 'ANTCGANT,ANTCGATC,ANTCCTNA,ANTCTTA,GATCGANT,GATCGATC,GATCCTNA,GATCTTA,TNAGGANT,TNAGGATC,TNAGCTNA,TNAGTTA,TAAGANT,TAAGATC,TAACTNA,TAATTA' \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney \
    -resume curious_pare
############3
# use this for background while out
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_R{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --digestion 'dpnii' \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney \
    -resume cranky_cajal


~/mnt/HiC-Pro-master/bin/utils/digest_genome.py -r G^ANTC ^GATC C^TNAG T^TAA -o pig_pg.bed /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

~/mnt/HiC-Pro-master/bin/HiC-Pro -i ~/mnt/muscle/ -o ~/mnt/eGen_HiC/phasegenomics/out/ -c ~/mnt/configMboI.txt 


/home/shiyi/mnt/hic/bin/mapped_2hic_fragments.py \
-f /home/shiyi/mnt/eGen_HiC/phasegenomics/pig_pg.bed \
-r /home/shiyi/mnt/hic/PG-WT-kidney/hicpro/mapping/sa22224_S3HiC_R.1_bwt2pairs.bam --all 


##############################
# downsampling
seqtk sample -s 123 /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_R1.fastq.gz 100000 | gzip > /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short1.fastq.gz
seqtk sample -s 123 /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_R2.fastq.gz 100000 | gzip > /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short2.fastq.gz

nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --dnase --min_cis 1000 \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney_downsample_dnase

nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --digestion 'arima' \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney_downsample_arima

nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --digestion 'dpnii' \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney_downsample_dpnii

nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney_downsample_4enzymes \
    --restriction_site 'G^ANTC ^GATC C^TNAG T^TAA' \
    --restriction_fragments /home/shiyi/mnt/eGen_HiC/phasegenomics/pig_pg.bed \
    --ligation_site 'ANTCGANT,ANTCGATC,ANTCCTNA,ANTCTTA,GATCGANT,GATCGATC,GATCCTNA,GATCTTA,TNAGGANT,TNAGGATC,TNAGCTNA,TNAGTTA,TAAGANT,TAAGATC,TAACTNA,TAATTA' 

#####
# adjusted ligatin product
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney_downsample_4enzymes_may18 \
    --restriction_site 'G^ANTC,^GATC,C^TNAG,T^TAA' \
    --restriction_fragments /home/shiyi/mnt/eGen_HiC/phasegenomics/pig_pg.bed \
    --ligation_site 'GANTANTC,GANTGATC,GANTTNAG,GANTTTAA,GATCANTC,GATCGATC,GATCTNAG,GATCTTAA,CTNAANTC,CTNAGATC,CTNATNAG,CTNATTAA,TTAAANTC,TTAAGATC,TTAATNAG,TTAATTAA' 
    
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-WT-kidney_downsample_4enzymes_may18 \
    --restriction_site 'G^ANTC ^GATC C^TNAG' \
    --ligation_site 'GANTANTC,GANTGATC,GANTTNAG,GANTTTAA,GATCANTC,GATCGATC,GATCTNAG,GATCTTAA,CTNAANTC,CTNAGATC,CTNATNAG,CTNATTAA,TTAAANTC,TTAAGATC,TTAATNAG,TTAATTAA' 
    
##################
# https://hicexplorer.readthedocs.io/en/latest/content/example_usage.html#reads-mapping

# multiple enzymes
# https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindRestSite.html#hicfindrestsite

################
# install
conda install hicexplorer -c bioconda -c conda-forge

################3
# align the reads

seqtk sample -s 123 /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_R1.fastq.gz 1000 | gzip > /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short1.fastq.gz
seqtk sample -s 123 /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_R2.fastq.gz 1000 | gzip > /home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short2.fastq.gz


bwa mem -A1 -B4  -E50 -L0  /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    sa22224_S3HiC_short1.fastq.gz 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam

bwa mem -A1 -B4  -E50 -L0  /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa     sa22224_S3HiC_short1.fastq.gz 2>>mate_R1.log >mate_R1.sam
bwa mem -A1 -B4  -E50 -L0  /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa     sa22224_S3HiC_short2.fastq.gz 2>>mate_R1.log >mate_R2.sam

bwa mem -A1 -B4  -E50 -L0  /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    sa22224_S3HiC_short2.fastq.gz 2>>mate_R2.log | samtools view -Shb - > mate_R2.bam

bwa mem -A1 -B4  -E50 -L0  /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa     sa22224_S3HiC_short2.fastq.gz 2>>mate_R1.log >mate_R2.sam
mate_R1.bam
mate_R1.bam

bwa mem -t 12 -A1 -B4  -E50 -L0  /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    sa22224_S3HiC_R1.fastq.gz 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam
bwa mem -t 12 -A1 -B4  -E50 -L0  /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    sa22224_S3HiC_R2.fastq.gz 2>>mate_R2.log | samtools view -Shb - > mate_R2.bam
###################
# cut the genome

hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern GAATC -o GAATC.bed
hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern GAGTC -o GAGTC.bed
hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern GATTC -o GATTC.bed
hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern GACTC -o GACTC.bed

hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern GATC -o GATC.bed
hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern CTAAG -o CTAAG.bed
hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern CTGAG -o CTGAG.bed
hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern CTTAG -o CTTAG.bed
hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern CTCAG -o CTCAG.bed

hicFindRestSite --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --searchPattern TTAA -o TTAA.bed

###############3
# bam to hic
# 2
hicBuildMatrix --samFiles mate_R1.sam mate_R2.sam \
                 --binSize 10000 \
                 --restrictionSequence GAATC GAGTC GATTC GACTC GATC CTAAG CTGAG CTTAG CTCAG TTAA \
                 --danglingSequence AAT AGT ATT ACT GATC TAA TGA TTA TCA TA \
                 --restrictionCutFile GAATC.bed GAGTC.bed GATTC.bed GACTC.bed GATC.bed CTAAG.bed CTGAG.bed CTTAG.bed CTCAG.bed TTAA.bed \
                 --threads 12 \
                 --inputBufferSize 100000 \
                 --outBam hic_2.bam \
                 -o hic_matrix_2.h5 \
                 -ga /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 --QCfolder ./hicQC_2

# testing using only GATC
# 3
hicBuildMatrix --samFiles mate_R1.bam mate_R2.bam \
                 --binSize 10000 \
                 --restrictionSequence GATC \
                 --danglingSequence GATC \
                 --restrictionCutFile GATC.bed \
                 --threads 10 \
                 --inputBufferSize 100000 \
                 --outBam hic_3.bam \
                 -o hic_matrix_3.h5 \
                 -ga /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 --QCfolder ./hicQC_3
# 1
hicBuildMatrix --samFiles mate_R1.sam mate_R2.sam \
                 --binSize 10000 \
                 --restrictionSequence GATC \
                 --danglingSequence GATC \
                 --restrictionCutFile GATC.bed \
                 --threads 12 \
                 --inputBufferSize 100000 \
                 --outBam hic.bam \
                 -o hic_matrix.h5 \
                 -ga /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 --QCfolder ./hicQC

#######################
hicCorrectMatrix diagnostic_plot -m hic_matrix.h5 -o hic_corrected.png






#########################
#################3
#############
################################
# reformat the dataframe
cut -f1,2,3,4,5,6,7,8,9,10 /home/shiyi/mnt/hic/PG-WT-kidney/hicpro/valid_pairs/sa22224_S3HiC_R.allValidPairs | \
awk '{OFS = "\t";print $1,$4,$2,$3,$9,$7,$5,$6,$10,$8}' | \
awk '{ OFS = "\t"; if($2=="+") $2=0 ; if($2=="-") $2=16 ;if($6=="+") $6=0 ; if($6=="-") $6=16 ; $1=""; $5=0; $9=1;print $0}' |  \
sort -k2,2d -k3,3n  | \
awk '{OFS = "\t";print $1,$2,$3,$4,$5,$6,$7,$8}'  > \
/home/shiyi/mnt/hic/PG_reformat.reorder.correct.sorted.pairs


#################
# making the hic file
# https://github.com/aidenlab/juicer/wiki/Pre
java -Xmx2g -jar /home/shiyi/mnt/myJuicerDir/scripts/common/juicer_tools.jar \
pre /home/shiyi/mnt/hic/PG_reformat.reorder.correct.sorted.pairs \
/home/shiyi/mnt/hic/PG_WT.hic \
/home/shiyi/mnt/2021genome/s11.chrom.sizes \
-d true

java -Xmx2g -jar ~/mnt/myJuicerDir/scripts/common/juicer_tools.jar hiccups \
--cpu --threads 14 \
-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000  \
--ignore_sparsity \
/home/shiyi/mnt/hic/PG_WT.hic \
/home/shiyi/mnt/hic/PG_WT.loops

###############################3
# using juicer
/home/shiyi/mnt/juicer/CPU/juicer.sh -h

Usage: juicer.sh [-g genomeID] [-d topDir] [-s site] [-a about] 
                 [-S stage] [-p chrom.sizes path] [-y restriction site file]
                 [-z reference genome file] [-D Juicer scripts directory]
                 [-b ligation] [-t threads] [-h] [-f] [-j]

# supply with ligation junction of dpnii
/home/shiyi/mnt/juicer/CPU/juicer.sh -d /home/shiyi/mnt/hic/PG-kidney_juicer_testing \
        -s DpnII -p /home/shiyi/mnt/2021genome/s11.chrom.sizes \
        -y /home/shiyi/mnt/hic/ss11_DpnII.txt \
                 -z /home/shiyi/mnt/2021genome/BWAIndex/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 -D /home/shiyi/mnt/eGen_HiC \
                 -t 12 

# no ligation junction of dpnii
/home/shiyi/mnt/juicer/CPU/juicer.sh -d /home/shiyi/mnt/hic/PG-kidney_juicer_testing \
        -p /home/shiyi/mnt/2021genome/s11.chrom.sizes \
        -y /home/shiyi/mnt/hic/ss11_DpnII.txt \
                 -z /home/shiyi/mnt/2021genome/BWAIndex/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 -D /home/shiyi/mnt/eGen_HiC \
                 -t 12 

# DpnII, Dde1, HinFI, and MseI
/home/shiyi/mnt/juicer/CPU/juicer.sh -d /home/shiyi/mnt/hic/PG-kidney_juicer_testing \
         -p /home/shiyi/mnt/2021genome/s11.chrom.sizes \
        -y /home/shiyi/mnt/hic/ss11_PG.txt \
                 -z /home/shiyi/mnt/2021genome/BWAIndex/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 -D /home/shiyi/mnt/eGen_HiC \
                 -t 12 
# this file format is wrong, i need to make new files

# no enzyme
/home/shiyi/mnt/juicer/CPU/juicer.sh -d /home/shiyi/mnt/hic/PG-kidney_juicer_testing \
         -p /home/shiyi/mnt/2021genome/s11.chrom.sizes \
         -s none \
                 -z /home/shiyi/mnt/2021genome/BWAIndex/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 -D /home/shiyi/mnt/eGen_HiC \
                 -t 12 
###########33
# real stuff
/home/shiyi/mnt/juicer/CPU/juicer.sh -d /home/shiyi/mnt/hic/PG-kidney_juicer \
         -p /home/shiyi/mnt/2021genome/s11.chrom.sizes \
        -y /home/shiyi/mnt/hic/ss11_PG.txt \
                 -z /home/shiyi/mnt/2021genome/BWAIndex/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 -D /home/shiyi/mnt/eGen_HiC \
                 -t 12 

#############################
# juicer format of re cut genome
python /home/shiyi/mnt/juicer/misc/generate_site_positions.py \
    DpnII ss11
python /home/shiyi/mnt/juicer/misc/generate_site_positions.py \
    PG ss11


######################################
# hayley loop calling
mustache -f /home/shiyi/mnt/pig_muscle_full.hic -ch 6 -r 5kb -pt 0.01 -o pig_muscle_chr6_out.tsv
mustache -f /home/shiyi/mnt/hic/PG_WT.hic -ch 6 -r 5kb -pt 0.01 -o pig_kidney_hicpro_chr6_out.tsv
mustache -f /home/shiyi/mnt/hic/PG-kidney_juicer/aligned/inter_30.hic -ch 6 -r 5kb -pt 0.01 -o pig_kidney_juicer_chr6_out.tsv
wc -l *.tsv
#####################
# the bug is magically solved
(mustache) shiyi@ip-192-168-1-231:~/mnt$ python -m mustache 
Traceback (most recent call last):
  File "/home/shiyi/miniconda3/envs/mustache/lib/python3.8/runpy.py", line 193, in _run_module_as_main
    return _run_code(code, main_globals, None,
  File "/home/shiyi/miniconda3/envs/mustache/lib/python3.8/runpy.py", line 86, in _run_code
    exec(code, run_globals)
  File "/home/shiyi/miniconda3/envs/mustache/lib/python3.8/site-packages/mustache/__main__.py", line 1, in <module>
    from .mustache import main
  File "/home/shiyi/miniconda3/envs/mustache/lib/python3.8/site-packages/mustache/mustache.py", line 14, in <module>
    import hicstraw
ModuleNotFoundError: No module named 'hicstraw'

#################################3
# helping Feng with genome assembly
aws s3 cp s3://egenesis-sagar/azenta_pacbio_HiFi/analysis/20114_with_HiC/
aws s3 cp s3://egenesis-sagar/azenta_pacbio_HiFi/analysis/20114_with_HiC/20114_hifi_reads.asm.hic.hap1.p_ctg.fa ./
nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/sa22224_S3HiC_short{1,2}.fastq.gz' \
    --fasta '/home/shiyi/mnt/eGen_HiC/20114_hifi_reads.asm.hic.hap1.p_ctg.fa' \
    --bwt2_index /home/shiyi/mnt/hic/bowtie2_index \
    --digestion 'dpnii' \
    --max_memory 32.GB \
    --max_cpus 12 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir PG-Yuc-WT-kidney_downsample_dpnii



aws s3 cp /home/shiyi/mnt/hic/PG-Yuc-WT-kidney_downsample_dpnii/hicpro/mapping/sa22224_S3HiC_short.1_bwt2pairs.bam s3://egenesis-data-processed/yin/genomic_assembly/



###########################################
###########################################
# A9161

aws s3 ls egenesis-sandbox/illumina/220720_VH00523_100_AAAK5G5HV/Analysis/1/Data/fastq/
aws s3 cp egenesis-sandbox/illumina/220720_VH00523_100_AAAK5G5HV/Analysis/1/Data/fastq/Undetermined_S0_L001_R1_001.fastq.gz
aws s3 cp egenesis-sandbox/illumina/220720_VH00523_100_AAAK5G5HV/Analysis/1/Data/fastq/Undetermined_S0_L001_R2_001.fastq.gz

####################################
# backend
###################################
sudo file -s /dev/nvme3n1
sudo mkfs -t xfs /dev/nvme4n1
# get file system
# sudo mount -t xfs -o nouuid /dev/nvme3n1 ~/mnt2
sudo mount /dev/nvme3n1 ~/mnt2

sudo chmod -R +777 mnt3
################################

seqtk sample -s 123 /home/shiyi/mnt2/A9161_hic/Undetermined_S0_L001_R1_001.fastq.gz 1000 | gzip > /home/shiyi/mnt/eGen_HiC/phasegenomics/A9161/30527_S3HiC_short1.fastq.gz
seqtk sample -s 123 /home/shiyi/mnt2/A9161_hic/Undetermined_S0_L001_R2_001.fastq.gz 1000 | gzip > /home/shiyi/mnt/eGen_HiC/phasegenomics/A9161/30527_S3HiC_short2.fastq.gz

/home/shiyi/mnt/eGen_HiC/phasegenomics/A9161/30527_S3HiC_short1.fastq.gz


nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt2/A9161_hic/Undetermined_S0_L001_R{1,2}_001.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --digestion 'dpnii' \
    --max_memory 32.GB \
    --max_cpus 10 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir /home/shiyi/mnt2/PG-A9161-kidney_downsample_dpnii \
    -resume kickass_watson \
    -bg


./HiC-Pro-master/bin/HiC-Pro -i ~/mnt/muscle/ -o ~/mnt/pigmuscleoutput/ -c ~/mnt/configMboI.txt 

cat /home/shiyi/mnt2/A9161_hic/Undetermined_S0_L00*_R1_001.fastq.gz > A9161_R1.fastq.gz
cat /home/shiyi/mnt2/A9161_hic/Undetermined_S0_L00*_R2_001.fastq.gz > A9161_R2.fastq.gz


nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt2/A9161_R{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --digestion 'dpnii' \
    --max_memory 32.GB \
    --max_cpus 50 \
    --split_fastq --fastq_chunks_size 200000 \
    --save_reference \
    --outdir /home/shiyi/mnt2/full_PG-A9161-kidney_downsample_dpnii \
#    -resume soggy_coulomb \
    -bg


nextflow run main.nf \
    -profile docker    \
    --input '/home/shiyi/mnt/eGen_HiC/phasegenomics/A9161/30527_S3HiC_short{1,2}.fastq.gz' \
    --bwt2_index /home/shiyi/mnt/2021genome/bowtie/   \
    --fasta /home/shiyi/mnt/2021genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
    --digestion 'dpnii' \
    --max_memory 32.GB \
    --max_cpus 10 \
    --split_fastq --fastq_chunks_size 2000000 \
    --save_reference \
    --outdir /home/shiyi/mnt2/testing_PG-A9161-kidney_downsample_dpnii \
    -bg

    /home/shiyi/mnt/eGen_HiC/phasegenomics/A9161/30527_S3HiC_short1.fastq.gz


docker pull nservant/hicpro:latest

docker login images.sbgenomics.com -u shiyi.yin@egenesisbio.com -p 64565d71daa343fbb9a9c919c4b4575f


The full run is failed

I can use manual deduplicated by picard
