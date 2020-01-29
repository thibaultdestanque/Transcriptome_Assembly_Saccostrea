#!/usr/bin/env nextflow


def helpMessage() {
  log.info"""

    =============================
      RNA-seq pipeline nextflow 
    =============================
    
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf  

    Tree folders should look like this: 
      DATA_TO_ANALYSE_FOLDER                  # Folder
          '=> main.nf                         # Script nextflow
          '=> nextflow.config                 # config nextflow
          '=> Ref_Genome                      # Folder
              '=> genome.fa                   # Reference genome in fasta
              '=> genomeIndex                 # Folder which will store index reference genome
                  '=> genome index files      # Files
              '=> genome_annotation.gff3      # Annotation of genome in gff3 
          '=> Fastq                           # Folder
              '=> fastq files to analyse      # Files to analysis in fastq
          '=> Other_files                     # Folder
              '=> adapter file                # adapter file

    Check presence of: 
     - genome.fa
     - genome_annotation.gff3
     - fastq files
     - adapter files
     - join_multiple_files script

  """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
 
// Show help message
params.help = false
if (params.help){
  helpMessage()
  exit 0
}


/*
 * Start RNA seq pipeline nextflow
 */
log.info ""
log.info ""
log.info "        ================================================="
log.info "              MAKE TRANSCRIPTOME PIPELINE NEXTFLOW "
log.info "        ================================================="
log.info ""
log.info ""
log.info "        Parameters used:"
log.info "            - fastq             : ${params.reads}"
log.info "            - Genome path       : ${params.genome_path}"
//log.info "            - Genome Annotation : ${params.annot}"
log.info ""
log.info "        Results will be store in:"
log.info "            ${params.outdir}"
log.info ""
log.info "        Process are stored in:"
log.info "            work/"
log.info ""
log.info ""
log.info "        ==================="
log.info "          Launch pipeline  "
log.info "        ==================="
log.info ""


/*
 * The reference genome file
 */
// genome_fa_file  = file(params.genome_fa)
// annotation_file = file(params.annot)
adapter_file    = file(params.adapter_file)
reads           = file(params.reads)



/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs } 

/*
 * Step 0. Pre-process fastqc
 */

 process Fastqc_row_data {
    publishDir "results/fastqc_row_data", mode: 'copy'

    time'2h'
    cpus 1
    queue 'sequentiel'
    memory '30 GB'
    echo true
    scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::fastqc'

    input:
    file(reads) from reads
     
    output:
    file "*" into fastqc

    shell:
    """
    fastqc ${reads} >& /home1/scratch/tdestanq/98_log_files/fastqc.log 2>&1
    ls $PWD/Fastq/*fastq.gz | awk '{print "echo \$(zcat "\$0" | wc -l)/4"}' | sh >& SeqCountFastq.txt 2>&1
    """
}
 
 /*
  * Step 1. Trim sequence with trimmomatic
  */
 
 process Trim {
     tag "$pair_id"
     publishDir "results/trimming", mode: 'copy'
 
     time'2h'
     cpus 8
     queue 'omp'
     memory '60 GB'
     echo true
     scratch '/home1/scratch/tdestanq/'
     conda 'bioconda::trimmomatic=0.36'
 
     input:
     set pair_id, file(reads) from read_pairs
     file(adapter_file) from adapter_file
      
     output:
     //file "*_paired.fastq.gz" into paired_read_trimed
     file "*_trim_R1_paired.fastq.gz" into paired_read_trimed_R1, paired_read_trimed_R1_for_trinity 
     file "*_trim_R2_paired.fastq.gz" into paired_read_trimed_R2, paired_read_trimed_R2_for_trinity 
     file "*_unpaired.fastq.gz" into unpaired_read_trimed
 
     shell:
     """
     trimmomatic PE -threads 8 -phred33 ${reads} \
     	${pair_id}_trim_R1_paired.fastq.gz \
     	${pair_id}_trim_R1_unpaired.fastq.gz \
     	${pair_id}_trim_R2_paired.fastq.gz \
     	${pair_id}_trim_R2_unpaired.fastq.gz \
     	ILLUMINACLIP:${adapter_file}:${params.illuminaclip_1}:${params.illuminaclip_2}:${params.illuminaclip_3} \
     	LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.slidingwindows_1}:${params.slidingwindows_2} MINLEN:${params.minlen} \
         >& /home1/scratch/tdestanq/98_log_files/Trim.log 2>&1
     """
 }

/*
 * Step 3. Trinity assembly
 */

 process Fastqc_trimmed {
    publishDir "results/fastqc_trimmed_data", mode: 'copy'

    time'2h'
    cpus 1
    queue 'sequentiel'
    memory '30 GB'
    echo true
    scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::fastqc'

    input:
    file(reads_R1) from paired_read_trimed_R1
    file(reads_R2) from paired_read_trimed_R2
     
    output:
    file "*" into fastqc_trimmed

    shell:
    """
    fastqc ${reads_R1} >& /home1/scratch/tdestanq/98_log_files/fastqc_trimmed_R1.log 2>&1
    fastqc ${reads_R2} >& /home1/scratch/tdestanq/98_log_files/fastqc_trimmed_R2.log 2>&1
    ls $PWD/results/trimming/*fastq.gz | awk '{print "echo \$(zcat "\$0" | wc -l)/4"}' | sh >& SeqCountFastq.txt 2>&1
    """

}


/*
 * Step 3. Trinity assembly
 */

process Trinity_Assembly {
    //publishDir "results/Trinity_Assembly", mode: 'copy'

    time'180h' // change it to 180h
    cpus 12
    queue 'omp'
    memory '160 GB'
    echo true
    //scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::trinity=2.8.5 bioconda::bowtie2 bioconda::salmon bioconda::samtools anaconda::numpy'
    
    input:
    file reads_R1 from paired_read_trimed_R1_for_trinity.toList()
    file reads_R2 from paired_read_trimed_R2_for_trinity.toList()

    output:
    file "*" into assembly_Trinity
    file "trinity_out_dir/Trinity.fasta" into trinity_fasta_for_BUSCO, trinity_fasta_for_Transdecoder


    shell:
    """
    bash
    ls *_trim_R1_paired.fastq.gz | tr '\n' ',' | awk 'sub(".\$","")' > LEFT.txt; # Add all R1 in a list separated by comas
    LEFT=`cat LEFT.txt`;
    ls *_trim_R2_paired.fastq.gz | tr '\n' ',' | awk 'sub(".\$","")' > RIGHT.txt;
    RIGHT=`cat RIGHT.txt`;

    conda list >& /home1/scratch/tdestanq/98_log_files/Env_list.log 2>&1 ;
    Trinity --seqType fq --max_memory 150G --left \$LEFT --right \$RIGHT --CPU 12 --verbose  >& /home1/scratch/tdestanq/98_log_files/Trinity.log 2>&1 ;
    """
}


/*
 * Step 4. BUSCO (Evaluate assembly) / Trinity stats
 */

process BUSCO_evaluate_assembly {
    publishDir "results/BUSCO_evaluate_assembly", mode: 'copy'

    time'2h' 
    cpus 1
    queue 'sequentiel'
    memory '50 GB'
    echo true
    //scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::busco=4.0.2 bioconda::trinity=2.8.5'

    
    input:
    file trinity_fasta_for_BUSCO from trinity_fasta_for_BUSCO

    output:
    file "*" into BUSCO_result

    shell:
    """
    busco -h >& busco_test.log 2>&1 ;
    head ${trinity_fasta_for_BUSCO} >& head_TrinityFa.txt 2>&1 ;
    busco -m transcriptome -i ${trinity_fasta_for_BUSCO} -o BUSCO_Transcriptome_Saccostrea -l /home1/datawork/tdestanq/BUSCO_downloads/eukaryota_odb10 >& BUSCOv4.0.1.log 2>&1 ;
    #TrinityStats.pl ${trinity_fasta_for_BUSCO} >& Trinity_stats.txt 2> /home1/scratch/tdestanq/98_log_files/Trinity_Stats.log
    """
}

/*
 * Step 5. Transdecoder (Predict ORF)
 */
 
process Transdecoder {
    publishDir "results/Transdecoder", mode: 'copy'

    time'5h' 
    cpus 1
    queue 'sequentiel'
    memory '50 GB'
    echo true
    scratch '/home1/scratch/tdestanq/'
    conda 'bioconda::transdecoder=5.5.0'
    
    input:
    file trinity_fasta_for_Transdecoder from trinity_fasta_for_Transdecoder

    output:
    file "*" into Transdecoder_results

    shell:
    """
    TransDecoder.LongOrfs -t ${trinity_fasta_for_Transdecoder} >& /home1/scratch/tdestanq/98_log_files/Transdecoder.longOrfs.log 2>&1 ;
    TransDecoder.Predict -t ${trinity_fasta_for_Transdecoder} >& /home1/scratch/tdestanq/98_log_files/Transdecoder.Predict.log 2>&1 ;
    """
}

/*
 * Step 6. TRINOTATE (Annotation)
 */

// process TRINOTATE_annotation {
//    publishDir "results/TRINOTATE_annotation", mode: 'copy'
//
//    time'5h' // change it to 180h
//    cpus 1
//    queue 'sequentiel'
//    memory '50 GB'
//    echo true
//    scratch '/home1/scratch/tdestanq/'
//    conda 'bioconda::trinotate=3.2.0'
//    
//    input:
//    file transdecoder_results from Transdecoder_results
//
//    output:
//    file "*" into TRINOTATE_annotation
//
//    shell:
//    """
//    
//    """
//}










workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Results files are stocked in --> $params.outdir ;)\n" : "Oops... something went wrong" )
}
