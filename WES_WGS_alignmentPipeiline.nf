/*
 * pipeline input parameters
*/

//projectDir = ""

params.reads = '/data/DaveJames/nextflow/RNASeqTutorial/training/nf-training/data/ggal/*_{1,2}.fq'
params.outdir = "/data/DaveJames/nextflow/WGS_pipeline/results/"
params.adaptorfile = "/data/genomes/GCF_000001405.40_BWA_MEM2_Index/adaptors.fa"
params.genomefile = "/data/genomes/G38_BROAD_150823/"
params.genomeid = "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

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
    
    publishDir params.outdir, mode: 'copy' 
    container 'rbgo/fastqc:0.12.0'

    input:
    tuple val(sample_id), path(reads)
    path results


    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    
    """


}


process SCYTHE {
    cpus 4
    container 'sickle_scythe_dj'

    input:
    path adaptor
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_Scythe{1,2}.fq.gz")
    
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

    gzip  ${sample_id}_Scythe2.fq
    gzip ${sample_id}_Scythe1.fq

    """


}

process SICKLE {
    cpus 4
    container 'sickle_scythe_dj'
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_ScytheSickle{1,2}.fq.gz")

    script:
    """
     sickle pe \
    -f ${reads[0]} \
    -r ${reads[1]} \
    -o ${sample_id}_ScytheSickle1.fq \
    -p  ${sample_id}_ScytheSickle2.fq \
    -s ${sample_id}_ScytheSickle_LowQCreads.fq -t sanger 

    gzip ${sample_id}_ScytheSickle1.fq
    gzip ${sample_id}_ScytheSickle2.fq
    gzip ${sample_id}_ScytheSickle_LowQCreads.fq
    """


}

process BWA_MEM2 {
    cpus 8 

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
    cpus 9
    //memory '20 GB'
    container 'biocontainers/samtools:v1.7.0_cv3'
    input:
    tuple val(sample_id), path(reads) 

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")
    //tuple val(sample_id), path("$")
    script:
    """
    samtools view -S -b ${reads[0]} > ${sample_id}.bam 
    samtools sort -@ 4 ${sample_id}.bam > ${sample_id}_sorted.bam 

    """

}


process MARK_DUPLICATES {
    cpus 10
    publishDir params.outdir, mode:'move'

    container 'broadinstitute/gatk:latest' 
    input:
    tuple val(sample_id), path(reads)
    output:
    path "${sample_id}_MarkedDup.bam" 
    path "${sample_id}_MarkedDuplicates.txt"
    path "${sample_id}_MarkedDup.bai"
    script:
    """
    gatk MarkDuplicates I=${reads[0]} O=${sample_id}_MarkedDup.bam M=${sample_id}_MarkedDuplicates.txt CREATE_INDEX=true 

    """   

}

process INDEXBAM {
    
    publishDir params.outdir, mode:'move'
    container 'biocontainers/samtools:v1.7.0_cv3'
    input:
    tuple val(sample_id), path(reads)
    path outDir    
    path bam
    path txt
    //output:
    
    //path "${sample_id}_MarkedDup.bam.bai"


    script:
    """
    
    samtools index -@ 12 ${outDir}/${sample_id}_MarkedDup.bam

    """

}

process ADDREADGROUPS {

    container 'broadinstitute/gatk:latest' 
    input:
    tuple val(sample_id), path(reads)

    output:
    
    tuple val(sample_id), path("${sample_id}_RG.bam")

    script:
    """

    gatk AddOrReplaceReadGroups I=${reads} O=${sample_id}_RG.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20


    """

 
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .view()
         .set { read_pairs_ch }
    
    fastqc_ch = FASTQC(read_pairs_ch, params.outdir)    
    scythe_ch = SCYTHE(params.adaptorfile, read_pairs_ch)
    //scythe_ch.view()
    sickle_ch = SICKLE(scythe_ch)
    //sickle_ch.view()
    bwa_ch = BWA_MEM2(sickle_ch, params.genomefile, params.genomeid)
    //bwa_ch.view()
    samtools_ch=SAM2BAM(bwa_ch)
    addreadgroups_ch=ADDREADGROUPS(samtools_ch)
    //samtools_ch.view()
    picard_ch=MARK_DUPLICATES(addreadgroups_ch)
//    picard_ch.view()
//    indexbam_ch=INDEXBAM(samtools_ch,params.outdir,picard_ch)
}
