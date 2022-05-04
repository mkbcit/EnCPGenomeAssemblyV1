=======================================================================================

       This bash script design for reference based chloroplast genome assembly
       written by Dr. Manosh Kumar Biswas @UOL  26April 2022  
 
       WorkFlow 1. Download Ensete NGS raw read from SRA database of  NCBI using SRA id   [Make SRA data list in text file]
                2. Check Quality then Clean and trim SRA reads
                3. Map clean R1 and R2 SRA reads in references CP genome [eg. Banana or ensete ]
                4. Extract consesus sequences assembly 
                5. Functional Annotation CP genes 
               
===============================================================================

   Dependent tools list: 1. sratoolkit [fastq-dump.3.0.0]
                         2. fastqc, multiQC, Trimmomatic 
                         3. BOWTIE2, samtools [view, bam2fq] 
                         
=============================================================================== 
