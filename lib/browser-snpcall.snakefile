shell.executable("/bin/bash")
# "unofficial bash strict mode" http://www.redsymbol.net/articles/unofficial-bash-strict-mode/
shell.prefix("set -euo pipefail;")

assert 'ACC' in config, ("required --config ACC=xyz missing")

# FIXME skipped all sed based source replacements in file
# FIXME add indel alignment quals and BAQ on the fly?

onerror:
    print("WORKFLOW ERROR")
    
onsuccess:
    print("WORKFLOW SUCCESS")


# first rule defines final target
rule final:
    input:
        gcov=config['ACC'] + ".gcov.gz",
        lofreq=config['ACC'] + ".lofreq.gz",
        srst2=config['ACC'] + ".srst2.gz",
#        lacer=config['ACC'] + ".lacer",
#        lacerdump=config['ACC'] + ".lacer.dump.gz"
    shell:
        "{config[DB_POST]} -mlst_database {config[SPECIES]} -source FASTQ {config[ACC]} -db;"


rule noref_final:
    input:
        tgz=config['ACC'] + ".tgz",
        srst2=config['ACC'] + ".srst2.gz"
    shell:
        "{config[DB_POST]} -mlst_database {config[SPECIES]} -source FASTQ {config[ACC]} -db;"


rule treefiles:
    input:
        gcov=config['ACC'] + ".gcov.gz",
        lofreq=config['ACC'] + ".lofreq.gz"


rule wgs:
    input:
        tgz=config['ACC'] + ".tgz"


rule cleanup:
    shell:
        "if {config[GETFILES]} -checkonly {config[ACC]}; then rm -f `{config[GETFILES]} {config[ACC]} -delimiter \" \"`; fi;"


rule bamonly:
    input:
        vsidq=config['ACC'] + "-vsidq.bam"
    output:
        bam=config['ACC'] + "-vsidq-full.bam"
    shell:
        "mv {input.vsidq} {output.bam} && {config[LOFREQ]} index {output.bam};"


rule lacer:
    input:
        ref=config['REF'],
        bam="{prefix}-vsidq.bam",
        bai="{prefix}-vsidq.bam.bai"
    output:
        lacer="{prefix}.lacer",
        lacerdump="{prefix}.lacer.dump.gz"
    threads:
        4
    shell:
        "{config[LACER]} -bam {input.bam} -ref {input.ref} -threads {threads} -outfile {output.lacer} -verbose 0 -gatk -dump 1 -stopbases {config[LACER_STOP]} | gzip > {output.lacerdump};"


rule srst2:
    output:
        srst2_file=config['ACC'] + ".srst2.gz",
        vffasta=config['ACC'] + ".resistance.fasta"
    threads:
        2
    shell:
        "set +u; "
        "source activate srst2; "
        "{config[SRST2WRAP]} -name {config[ACC]} -q1 `{config[GETFILES]} {config[ACC]} -delimiter \" -q2 \"` -species {config[SPECIES]} | gzip > {output.srst2_file}; "
        "source deactivate; "
        "set -u; "
        "{config[VFFASTA]} {output.srst2_file} > {output.vffasta};"


rule map:
    input:
        ref=config['REF']
    output:
        bam=temp(config['ACC'] + ".bam")
    threads:
        8
    shell:
        "if [ `{config[GETFILES]} {config[ACC]} | wc -l` -ge 1 ]; then "
        "  {config[BWA]} mem -t {threads} {input.ref} `{config[GETFILES]} {config[ACC]} -delimiter \" \"` | {config[SAMTOOLS]} view -@ {threads} -bS - > {output.bam}; "
        "else "
        "  false; "
        "fi;"
        

rule viterbi:
    input:
        ref=config['REF'],
        bam="{acc}.bam"
    threads:
        8
    output:
        bam=temp("{acc}-vsidq.bam")
    shell:
        "{config[LOFREQ]} viterbi -f {input.ref} {input.bam} | {config[LOFREQ]} indelqual --dindel -f {input.ref} /dev/stdin | {config[LOFREQ]} alnqual -b -r - {input.ref} | {config[SAMTOOLS]} sort -T {config[ACC]} -@ {threads} -O bam -o {output.bam} -;"

        
rule bam_index:
    input:
        bam="{prefix}.bam"
    output:
        bam=temp("{prefix}.bam.bai")
    shell:
        "{config[LOFREQ]} index {input.bam}"


rule lofreq_call:
    input:
        ref = config["REF"],
        bam = "{prefix}-vsidq.bam",
        bai = "{prefix}-vsidq.bam.bai"
    threads:
        8
    output:
        vcf = temp("{prefix}.lofreq")
    shell:
        "{config[LOFREQ]} call-parallel --pp-threads {threads} --call-indels -f {input.ref} -o {output.vcf} {input.bam};"
        
        
rule vcf_clean:
    input:
        lf_vcf="{acc}.lofreq"
    output:
        lf_vcf_clean="{acc}.lofreq.gz",
        lf_vcf_clean_tbi="{acc}.lofreq.gz.tbi"
    message: "Cleaning {input.lf_vcf} to create {output.lf_vcf_clean}"
    shell:
        "{config[LFVCF]} {input.lf_vcf} -sample {config[ACC]} | {config[BGZIP]} -c > {output.lf_vcf_clean} && {config[TABIX]} -p vcf {output.lf_vcf_clean}"


rule gcov:
    input:
        bam="{acc}-vsidq.bam",
        ref=config['REF']
    output:
        gcovgz="{acc}.gcov.gz",
        gcovpng="{acc}.gcov.png"
    params:
        # used throughout file to append. can therefore not be input or output of any rule
        gcov = lambda wildcards: wildcards.acc + ".gcov"
    message: "Running gcov on {input.bam} to generate {output.gcovgz}"
    shell:
        "if file --mime-type `{config[GETFILES]} {config[ACC]} | head -n 1` | grep -q \"gzip\$\"; then "
        "  zcat `{config[GETFILES]} {config[ACC]} | head -n 1` | wc -l | sed -e \"s/^/# {config[ACC]} fastq lines: /\" >> {params.gcov}; "
        "else "
        "  cat `{config[GETFILES]} {config[ACC]} | head -n 1` | wc -l | sed -e \"s/^/# {config[ACC]} fastq lines: /\" >> {params.gcov}; "
        "fi; "
        "{config[SAMTOOLS]} view -H {input.bam} | sed -e 's/^/# /' >> {params.gcov}; "
        "{config[BEDTOOLS]} genomecov -ibam {input.bam} -d -g {input.ref} >> {params.gcov}; "
        "gzip -f {params.gcov}; "
        "{config[GCOVGRAPH]} {params.gcov}.gz;"


rule assembly:
    output:
        wgstgz=config['ACC'] + ".tgz",
        wgsgbk=config['ACC'] + ".gbk"
    threads:
        8
    shell:
        "{config[GERMS_WGS]} -q1 `{config[GETFILES]} {config[ACC]} -delimiter \" -q2 \"` -output_dir . -t {threads} -species \"{config[FULLSPECIES]}\" -name {config[ACC]} && tar xvzf {config[ACC]}.tgz --strip-components=1 {config[ACC]}/{config[ACC]}.gbk;"
