#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//
// WORKFLOW
//
workflow {

    reads_ch = Channel.fromFilePairs(params.fastq, flat: true)
    ref_ch   = Channel.value(file(params.reference))

    fastqc_out = FASTQC(reads_ch)
    trimmed    = TRIMMOMATIC(reads_ch)
    aligned    = BWA_MEM(trimmed, ref_ch)
    bam_file   = SAMTOOLS_SORT(aligned)
    vcf_final  = BCFTOOLS_CALL(bam_file, ref_ch)

    db_file      = VCF_TO_SQLITE(vcf_final)
    summary_file = VARIANT_SUMMARIES(db_file)

    println "Final VCF: ${params.outdir}/vcf/variants.vcf"
    println "SQLite DB: ${params.outdir}/db/variants.db"
    println "Summary CSV files:"
    println "  - results/db/variant_summary.csv"
    println "  - results/db/snv_indel_counts.csv"
    println "  - results/db/substitution_spectrum.csv"
    println "  - results/db/top_quality_variants.csv"
}

//////////////////////////////////////////////////////////////
// FASTQC
//////////////////////////////////////////////////////////////
process FASTQC {
    tag "FASTQC"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
        tuple val(sample), path(read1), path(read2)

    output:
        path "*.html"
        path "*.zip"

    script:
    """
    fastqc ${read1} ${read2}
    """
}

//////////////////////////////////////////////////////////////
// TRIMMOMATIC
//////////////////////////////////////////////////////////////
process TRIMMOMATIC {
    tag "TRIMMOMATIC"
    container "trimmomatic_user.sif"
    publishDir "${params.outdir}/trimmomatic", mode: 'copy'

    input:
        tuple val(sample), path(read1), path(read2)

    output:
        tuple path("R1_paired.fq.gz"),
              path("R2_paired.fq.gz")

    script:
    """
    java -jar \$TRIM PE \
        ${read1} ${read2} \
        R1_paired.fq.gz R1_unpaired.fq.gz \
        R2_paired.fq.gz R2_unpaired.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25
    """
}

//////////////////////////////////////////////////////////////
// BWA MEM
//////////////////////////////////////////////////////////////
process BWA_MEM {
    tag "BWA_MEM"
    publishDir "${params.outdir}/bwa", mode: 'copy'

    input:
        tuple path(r1), path(r2)
        path ref

    output:
        path "aligned.sam"

    script:
    """
    bwa index ${ref}
    bwa mem ${ref} ${r1} ${r2} > aligned.sam
    """
}

//////////////////////////////////////////////////////////////
// SAMTOOLS SORT
//////////////////////////////////////////////////////////////
process SAMTOOLS_SORT {
    tag "SAMTOOLS"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
        path sam

    output:
        path "sorted.bam"

    script:
    """
    samtools view -bS ${sam} | samtools sort -o sorted.bam
    """
}

//////////////////////////////////////////////////////////////
// BCFTOOLS CALL
//////////////////////////////////////////////////////////////
process BCFTOOLS_CALL {
    tag "BCFTOOLS"
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
        path bam
        path ref

    output:
        path "variants.vcf"

    script:
    """
    samtools index ${bam}
    bcftools mpileup -f ${ref} ${bam} | \
        bcftools call -mv -Ov -o variants.vcf
    """
}

//////////////////////////////////////////////////////////////
// VCF → SQLite 
//////////////////////////////////////////////////////////////
process VCF_TO_SQLITE {
    tag "VCF_TO_SQLITE"
    publishDir "${params.outdir}/db", mode: 'copy'

    input:
        path vcf

    output:
        path "variants.db"

    script:
    """
    python3 << 'EOF'
import sqlite3

db = sqlite3.connect("variants.db")
cur = db.cursor()

# Create table
cur.execute("DROP TABLE IF EXISTS variants")
cur.execute(
    "CREATE TABLE variants ("
    " chromosome TEXT,"
    " position INTEGER,"
    " reference TEXT,"
    " alternative TEXT,"
    " quality REAL)"
)

# Parse VCF
with open(r"${vcf}", "r") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\\n").split("\\t")
        if len(cols) < 6:
            continue
        chrom, pos, _id, ref, alt, qual = cols[:6]
        try:
            pos = int(pos)
        except:
            continue
        if qual == "." or qual == "":
            qval = None
        else:
            try:
                qval = float(qual)
            except:
                qval = None
        cur.execute(
            "INSERT INTO variants VALUES (?,?,?,?,?)",
            (chrom, pos, ref, alt, qval)
        )

db.commit()
db.close()
EOF
    """
}

//////////////////////////////////////////////////////////////
// Additional Results (CSV files from the SQLite DB)
//////////////////////////////////////////////////////////////
process VARIANT_SUMMARIES {
    tag "VARIANT_SUMMARIES"
    publishDir "${params.outdir}/db", mode: 'copy'

    input:
        path db

    output:
        path "variant_summary.csv"
        path "snv_indel_counts.csv"
        path "substitution_spectrum.csv"
        path "top_quality_variants.csv"

    script:
    """
    python3 << 'EOF'
import sqlite3, csv

con = sqlite3.connect(r"${db}")
cur = con.cursor()

# 1) Full list by position
rows = list(cur.execute(
    "SELECT chromosome, position, reference, alternative, quality "
    "FROM variants "
    "ORDER BY position"
))
with open("variant_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["chromosome","position","reference","alternative","quality"])
    w.writerows(rows)

# 2) SNVs vs INDELs
rows = list(cur.execute(
    "SELECT CASE WHEN length(reference)=1 AND length(alternative)=1 "
    "THEN 'SNV' ELSE 'INDEL' END AS type, "
    "COUNT(*) AS count "
    "FROM variants "
    "GROUP BY type "
    "ORDER BY count DESC"
))
with open("snv_indel_counts.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["type","count"])
    w.writerows(rows)

# 3) Substitution spectrum (SNVs only)
rows = list(cur.execute(
    "SELECT reference || '>' || alternative AS substitution, COUNT(*) AS count "
    "FROM variants "
    "WHERE length(reference)=1 AND length(alternative)=1 "
    "GROUP BY substitution "
    "ORDER BY count DESC"
))
with open("substitution_spectrum.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["substitution","count"])
    w.writerows(rows)

# 4) Top 10 by quality
rows = list(cur.execute(
    "SELECT chromosome, position, reference, alternative, quality "
    "FROM variants "
    "ORDER BY quality DESC "
    "LIMIT 10"
))
with open("top_quality_variants.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["chromosome","position","reference","alternative","quality"])
    w.writerows(rows)

con.close()
EOF
    """
}