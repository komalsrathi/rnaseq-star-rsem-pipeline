from util.varsub import varsub

# config file 
configfile: "config.yaml"
varsub(config)

shell.prefix("source ~/.bash_profile")

# determine which genome reference you would like to use
# here we are using mm10 + ERCC92
# depending on the freeze, the appropriate references and data files will be chosen from the config
freeze = config['freeze']

# read list of samples, one per line
with open(config['datadirs']['samples']) as f:
    SAMPLES = f.read().splitlines()

rule all:
    input:
        starindex = config['reference']['stargenomedir'][freeze] + "/" + "SAindex",
        rsemindex = config['reference']['rsemgenomedir'][freeze] + ".n2g.idx.fa",
        fastqs = expand(config['sourcedir'] + config['datadirs']['fastq'] + "/" + "{file}_{rep}.fastq.gz", file = SAMPLES, rep = ['1','2']),
        bams = expand(config['sourcedir'] + config['datadirs']['bam'] + "/" + "{file}_Aligned.toTranscriptome.out.bam", file = SAMPLES),
        genes = expand(config['sourcedir'] + config['datadirs']['quant'] + "/" + "{file}.genes.results", file = SAMPLES)

# convert SRA to fastq
# this will be run only if you are starting with .sra files 
rule fastq_dump:
    input:
        sra = config['sourcedir'] + config['datadirs']['sra'] + "/" + "{file}.sra"
    params:
        fastqdump = config['tools']['fastqdump'],
        outdir = config['sourcedir'] + config['datadirs']['fastq']
    threads: 10
    output:
        f1 = config['sourcedir'] + config['datadirs']['fastq'] + "/" + "{file}_1.fastq.gz",
        f2 = config['sourcedir'] + config['datadirs']['fastq'] + "/" + "{file}_2.fastq.gz"
    shell:
        """
        {params.fastqdump} --split-3 --gzip --outdir {params.outdir} {input.sra}
        """

# create STAR genome index
rule star_genome:
    input:
        fasta = config['reference']['fasta'][freeze],
        gtf = config['reference']['gtf'][freeze]
    output:
        genomedir = config['reference']['stargenomedir'][freeze],
        starindex = config['reference']['stargenomedir'][freeze] + "/" + "SAindex"
    params:
        star = config['tools']['star'],
        overhang = 99
    threads: 4
    shell:
        """
        mkdir -p {output.genomedir}

        {params.star} \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output.genomedir} \
        --outFileNamePrefix {output.genomedir}
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.overhang}
        """

# create RSEM genome index
rule rsem_genome:
    input:
        fasta = config['reference']['fasta'][freeze],
        gtf = config['reference']['gtf'][freeze],
        genomedir = config['reference']['rsemgenomedir'][freeze]
    output:
        rsemindex = config['reference']['rsemgenomedir'][freeze] + ".n2g.idx.fa"
    params:
        prepref = config['tools']['rsem']['prepref']
    threads: 4
    shell:
        """
        {params.prepref} \
        -p {threads} \
        --gtf {input.gtf} {input.fasta} {input.genomedir}
        """

# align using STAR
rule star_align:
    input:
        f1 = config['sourcedir'] + config['datadirs']['fastq'] + "/" + "{file}_1.fastq.gz",
        f2 = config['sourcedir'] + config['datadirs']['fastq'] + "/" + "{file}_2.fastq.gz"
    output:
        out = config['sourcedir'] + config['datadirs']['bam'] + "/" + "{file}_Aligned.toTranscriptome.out.bam"
    params:
        star = config['tools']['star'],
        genomedir = config['reference']['stargenomedir'][freeze],
        prefix = config['sourcedir'] + config['datadirs']['bam'] + "/" + "{file}_"
    threads: 8
    shell:  
        """
        {params.star} \
        --runThreadN {threads} \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.f1} {input.f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype None \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM \
        --outSAMattributes NH HI AS NM MD \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --limitBAMsortRAM 50000000000
        """

# quantify expression using RSEM
rule rsem_norm:
    input:
        bam = config['sourcedir'] + config['datadirs']['bam'] + "/" + "{file}_Aligned.toTranscriptome.out.bam"
    output:
        genes = config['sourcedir'] + config['datadirs']['quant'] + "/" + "{file}.genes.results"
    params:
        calcexp = config['tools']['rsem']['calcexp'],
        genomedir = config['reference']['rsemgenomedir'][freeze],
        prefix = config['sourcedir'] + config['datadirs']['quant'] + "/" + "{file}"
    threads: 6
    shell:
        """
        {params.calcexp} \
        --paired-end \
        --no-bam-output \
        --quiet \
        --no-qualities \
        -p {threads} \
        --forward-prob 0.5 \
        --seed-length 25 \
        --fragment-length-mean -1.0 \
        --bam {input.bam} {params.genomedir} {params.prefix}
        """

