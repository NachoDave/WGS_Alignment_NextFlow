/*
 * pipeline input parameters
*/

//projectDir = ""

params.reads = '/data/DaveJames/nextflow/RNASeqTutorial/training/nf-training/data/ggal/*_{1,2}.fq'
params.genome = ""

params.outdir = "/data/DaveJames/nextflow/WGS_pipeline/results/"
params.adaptorfile = "/data/genomes/GCF_000001405.40_BWA_MEM2_Index/adaptors.fa"
params.genomefile = "/data/genomes/GCF_000001405.40_BWA_MEM2_Index/"
params.genomeid = "GCF_000001405.40_GRCh38.p14_genomic.fna"
log.info """\

    RBGO WES/WGS -NF PIPELINE
    =========================
reads : ${params.reads}

"""
 .stripIndent()

/*
 * FASTQC process which runs fastqc on all of the fastq files
*/

process FASTQC {
    container 'rbgo/fastqc:0.12.0'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}

    """


}


process SCYTHE {

    container 'sickle_scythe_dj'

    input:
    path adaptor
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_Scythe{1,2}.fq")
    
    script:
    """
    scythe -a $adaptor \
    -o ${sample_id}_Scythe1.fq \
    -m ${sample_id}_matches1.txt \
    ${reads[0]} 

    scythe -a $adaptor \
    -o ${sample_id}_Scythe2.fq \
    -m ${sample_id}_matches2.txt \
    ${reads[1]}
    """

}

process SICKLE {

    container 'sickle_scythe_dj'
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_ScytheSickle{1,2}.fq")

    script:
    """
     sickle pe \
    -f ${reads[0]} \
    -r ${reads[1]} \
    -o ${sample_id}_ScytheSickle1.fq \
    -p  ${sample_id}_ScytheSickle2.fq \
    -s ${sample_id}_ScytheSickle_LowQCreads.fq -t sanger 

    """


}

process BWA_MEM2 {
    cpus 4

    container 'rbgo/bwa-mem2:2.2.1_v1'
    input:
    tuple val(sample_id), path(reads)
    path(genome)
    val genomeid

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    bwa-mem2 mem -t 4 ${genome}/$genomeid ${reads[0]} ${reads[1]} > ${sample_id}.sam

    """


}

process SAM2BAM {
    cpus 4

    container 'biocontainers/samtools:v1.7.0_cv3'
    input:
    tuple val(sample_id), path(reads) 

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    script:
    """
    samtools view -S -b ${reads[0]} | samtools sort -@ 10 -m 5G > ${sample_id}_sorted.bam 
    samtools index -@ 12 ${sample_id}_sorted.bam

    """

}


process MARK_DUPLICATES {
    cpus 4
    publishDir params.outdir, mode:'move'

    container 'broadinstitute/gatk:latest' 
    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_MarkedDup.bam"
    path "${sample_id}_MarkedDuplicates.txt"

    script:
    """
    gatk MarkDuplicates I=${reads[0]} O=${sample_id}_MarkedDup.bam M=${sample_id}_MarkedDuplicates.txt

    """   

}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .view()
         .set { read_pairs_ch }
    
    fastqc_ch = FASTQC(read_pairs_ch)    
    scythe_ch = SCYTHE(params.adaptorfile, read_pairs_ch)
    scythe_ch.view()
    sickle_ch = SICKLE(scythe_ch)
    sickle_ch.view()
    bwa_ch = BWA_MEM2(sickle_ch, params.genomefile, params.genomeid)
    bwa_ch.view()
    samtools_ch=SAM2BAM(bwa_ch)
    samtools_ch.view()
    picard_ch=MARK_DUPLICATES(samtools_ch)
    //picard_ch.view()
}
