/*
 * pipeline input parameters
*/

//projectDir = ""

params.reads = '/data/DaveJames/nextflow/RNASeqTutorial/training/nf-training/data/ggal/gut_{1,2}.fq'
params.genome = ""

params.outdir = "$projectDir/results/"
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

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .view()
         .set { read_pairs_ch }
    
//fastqc_ch = FASTQC(read_pairs_ch)    
    scythe_ch = SCYTHE(params.adaptorfile, read_pairs_ch)
    scythe_ch.view()
    sickle_ch = SICKLE(scythe_ch)
    sickle_ch.view()
    bwa_ch = BWA_MEM2(sickle_ch, params.genomefile, params.genomeid)
    bwa_ch.view()
    

}
