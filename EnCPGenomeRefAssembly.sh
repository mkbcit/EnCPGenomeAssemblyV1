#!/bin/sh

#=====================================================================================================================================
#       This bash script design for reference based chloroplast genome assembly
#       written by Dr. Manosh Kumar Biswas @UOL  26April 2022  
# 
#       WorkFlow 1. download SRA data from NCBI using SRA id   [SRA data list and details in text file]
#                2. Check Quality then Clean and trim SRA reads
#                3. Map clean R1 and R2 SRA reads in references CP genome [eg. Banana or ensete ]
#                4. Extract consesus sequences assembly 
#                5. Functional Annotation CP genes 
#               
#======================================================================================================================================
#####   Dependent tools list: 1. sratoolkit [fastq-dump.3.0.0]
#####                         2. fastqc, multiQC, Trimmomatic 
#####                         3. BOWTIE2, samtools [view, bam2fq] 
#=======================================================================================================================================  


##### Step 1 download data from NCBI R1 and R2 formate
    	mkdir RAW
  	cat EnSRAlist.txt | parallel "/home/mkb35/sratoolkit/bin/fastq-dump.3.0.0  --gzip --split-files {} -O RAW/" ###### zip formate download ###
  	
##### Step 2 Quality Control checking using fastQC quality and multiQuality 
 	mkdir RAWQC
 	fastqc -t 64 RAW/*.fastq.gz -O RAWQC/
	
	#### NOTE:   to run multiQC we need to active py3.7
 	conda activate py3.7  
	
         multiqc RAWQC/.

	
##### Step 3 Cleaning Triming Adapter seq,$INPUTFILE/ low quality reads using TrimmomaticPE
 	mkdir CleanRead
        cat EnSRAlist.txt | \
 	parallel "TrimmomaticPE     -threads  30    -phred33   RAW/{}_1.fastq.gz   RAW/{}_2.fastq.gz  \
 	CleanRead/clean{}_1.fastq.gz  CleanRead/U{}_1.fastq.gz   CleanRead/clean{}_2.fastq.gz    CleanRead/U{}_2.fastq.gz  \
 	-summary   CleanRead/TrimSummary{}.txt  \
 	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 \
 	TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"

##### Step 4 Index build for BOWTIE2 mapping pipeline
########  refgenome path: /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/ 
       bowtie2-build  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/RefEnvCP.fasta   /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome
       
##### Step5  Bowtie2 align         
	INPUTFILE=/media/mkb35/Madatacenter/EnseteNCBIRead/CleanRead
    
 	echo " mapping ........    Arkiya \n"
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                -1 $INPUTFILE/cleanSRR4304988_1.fastq.gz,$INPUTFILE/cleanSRR4304981_1.fastq.gz,$INPUTFILE/cleanSRR4304987_1.fastq.gz,$INPUTFILE/cleanSRR4304969_1.fastq.gz,$INPUTFILE/cleanSRR4304970_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304988_2.fastq.gz,$INPUTFILE/cleanSRR4304981_2.fastq.gz,$INPUTFILE/cleanSRR4304987_2.fastq.gz,$INPUTFILE/cleanSRR4304969_2.fastq.gz,$INPUTFILE/cleanSRR4304970_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Arkiya.sam		
         
	echo " mapping ........    Astara \n"
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1   $INPUTFILE/cleanSRR4304989_1.fastq.gz \
                 -2   $INPUTFILE/cleanSRR4304989_2.fastq.gz \
                 -S   /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Astara.sam 		
 	
	echo " mapping ........    Bedadit   \n"
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR1515268_1.fastq.gz,$INPUTFILE/cleanSRR1515269_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR1515268_2.fastq.gz,$INPUTFILE/cleanSRR1515269_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Bedadit.sam 		

       echo " mapping ........    Buffero  \n"
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304990_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304990_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Buffero.sam 		
       	

	echo " mapping ........    C Barrett 359 CA    \n" 
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR2923345_1.fastq.gz,$INPUTFILE/cleanSRR7439734_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR2923345_2.fastq.gz,$INPUTFILE/cleanSRR7439734_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/CBarrett359CA.sam 		

	echo " mapping ........      Derea   \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4308285_1.fastq.gz,$INPUTFILE/cleanSRR4308286_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4308285_2.fastq.gz,$INPUTFILE/cleanSRR4308286_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Derea.sam 		
 
	echo " mapping ........    Erpha_13  \n"
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304991_1.fastq.gz,$INPUTFILE/cleanSRR4304992_1.fastq.gz,$INPUTFILE/cleanSRR4304971_1.fastq.gz,$INPUTFILE/cleanSRR4304993_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304991_2.fastq.gz,$INPUTFILE/cleanSRR4304992_2.fastq.gz,$INPUTFILE/cleanSRR4304971_2.fastq.gz,$INPUTFILE/cleanSRR4304993_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Erpha.sam 		
    
	echo " mapping ........    Gena  \n"
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR974726_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR974726_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Gena.sam 

	echo " mapping ........    India  \n"
       bowtie2 --local  --threads 30 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR8879209_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR8879209_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/India.sam 
	
	echo " mapping ........    JungleSeeds   \n" 
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR610420_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR610420_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/JungleSeeds.sam 

	echo " mapping ........    Lochingie  \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304974_1.fastq.gz,$INPUTFILE/cleanSRR4304975_1.fastq.gz,$INPUTFILE/cleanSRR4304972_1.fastq.gz,$INPUTFILE/cleanSRR4304973_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304974_2.fastq.gz,$INPUTFILE/cleanSRR4304975_2.fastq.gz,$INPUTFILE/cleanSRR4304972_2.fastq.gz,$INPUTFILE/cleanSRR4304973_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Lochingie.sam 
	
	echo " mapping ........    Mazia  \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304976_1.fastq.gz,$INPUTFILE/cleanSRR4304977_1.fastq.gz,$INPUTFILE/cleanSRR4304978_1.fastq.gz,$INPUTFILE/cleanSRR4304979_1.fastq.gz,$INPUTFILE/cleanSRR4304980_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304976_2.fastq.gz,$INPUTFILE/cleanSRR4304977_2.fastq.gz,$INPUTFILE/cleanSRR4304978_2.fastq.gz,$INPUTFILE/cleanSRR4304979_2.fastq.gz,$INPUTFILE/cleanSRR4304980_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Mazia.sam 

	echo " mapping ........    Nechuwe     \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304982_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304982_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Nechuwe.sam 

	echo " mapping ........    Nobo  \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304983_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304983_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Nobo.sam 

	echo " mapping ........    Onjamo  \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4308284_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4308284_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Onjamo.sam 

	echo " mapping ........    Siyuti  \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304984_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304984_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Siyuti.sam 

	echo " mapping ........    Yako   \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304985_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304985_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Yako.sam 
     
	echo " mapping ........    Yanbule  \n"
       bowtie2 --local  --threads 35 -x  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/enseterefcpgenome \
                 -1 $INPUTFILE/cleanSRR4304986_1.fastq.gz \
                 -2 $INPUTFILE/cleanSRR4304986_2.fastq.gz \
                 -S /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Yanbule.sam 	
                 
        echo "\n---- Bowtie2 mapping finish----\n-----SAM conversion to BAM Start------\n"         

##### Step6  convert the SAM file into a BAM file and sort     
	mkdir CPconsesus_seq
	cat cpmapsamlist.txt | parallel "samtools view -bS /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/{}.sam | samtools sort - -o /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Sort{}.bam"
	
######## Step7  extract consesus sequences  in fastq fromate 
	cat cpmapsamlist.txt | parallel "samtools mpileup -uf /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/RefEnvCP.fasta  /media/mkb35/Madatacenter/EnseteNCBIRead/En15CP_3rdround_analysis/Sort{}.bam | bcftools call -c | vcfutils.pl vcf2fq > CPconsesus_seq/Consen_CP_{}.fastq"
 
####### Step8 Convert .fastq to .fasta and set bases of quality lower than 20 to N
	cat cpmapsamlist.txt | parallel "seqtk seq -aQ64 -q20 -n N CPconsesus_seq/Consen_CP_{}.fastq > CPconsesus_seq/Consen_CP_{}.fasta"              
              
####### Step9 rename id with file name
	for file in *.fasta;
	   do
	   sed -i "s/>.*/${file%%.*}/" "$file" ;
       done  
              
####### Step10 coppy all the seq in a single file
        cat *.fasta >> allcpConsesu.fasta
              

