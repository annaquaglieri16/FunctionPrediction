nextflow.enable.dsl = 2

process FETCH {
    publishDir "${params.outdir}/data/raw", mode: 'copy'

    input:
    path organism_config

    output:
    path "proteome.fasta"

    script:
    """
    plm-fetch --organism-config $organism_config --out-path proteome.fasta
    """
}

process DIGEST {
    publishDir "${params.outdir}/processed", mode: 'copy'

    input:
    path fasta
    path enzyme_config

    output:
    path "peptides.parquet"

    script:
    """
    plm-digest --fasta-path $fasta --enzyme-config $enzyme_config --out-path peptides.parquet
    """
}

workflow {
    organism_cfg = Channel.fromPath(params.organism_config)
    enzyme_cfg   = file(params.enzyme_config)

    FETCH(organism_cfg)
    DIGEST(FETCH.out, enzyme_cfg)
}
