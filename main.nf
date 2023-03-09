def helpMessage() {
        log.info"""
Usage:
Mandatory arguments:
Velvet options:
BUSCO options:
Run options:
Skip metrics:
    """.stripIndent()
}

params.index = file(params.input)
params.outdir = "results"
params.phylo_method = "fasttree"

workflow {
    Boolean flag = false
    Channel.fromPath(params.index) \
    .ifEmpty {exit 1, log.info "Cannot find path file ${csvFile}"} \
    .splitText {it.strip() } \
    .map { it -> tuple(it) } \
    .set {axt}
    main:
     split(axt)
     align(split.out)
     trim(align.out)
     concatx(trim.out.collect())
    if (params.phylo_method == "fasttree" ) {
    fasttree(concatx.out)
    }
    if (params.phylo_method == "iqtree" ) {
    iqtree(concatx.out)
    }
}

process split {
    cpus = 4
    tag "X"

    input:
    tuple val(genes)

    output:
    tuple val(genes), path("${genes}.fa")

    script:
    """
    seqkit grep -r -n -p ${genes} $baseDir/subset.test.fasta > ${genes}.fa
    """
}

process align {
    cpus = 1
    tag "X"

    input:
    tuple val(genes), path("${genes}.fa")

    output:
    tuple val(genes), path("${genes}.aln")

    script:
    """
    mafft --thread 1 --auto ${genes}.fa > ${genes}.aln
    """
}

process trim {
    cpus = 1
    tag "X"

    input:
    tuple val(genes), path("${genes}.aln")

    output:
    path("${genes}.trimal")

    script:
    """
    trimal -in ${genes}.aln -out ${genes}.trimal
    """
}

process concatx {
    cpus = 1
    tag "X"
    publishDir "$params.outdir/", mode: 'copy'

    input:
    path("*.trimal")

    output:
    path("concat_all.txt")

    script:
    """
    python $baseDir/scripts/concat_aln.py --input *.trimal --output concat_all.txt
    """
}

process fasttree {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("concat_all.txt")

    output:
    path("out.tree")

    script:
    """
    FastTreeMP -gtr -nt -log logfile < concat_all.txt > out.tree
    """
}

process iqtree {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("concat_all.txt")

    output:
    path("concat_all.txt.treefile")

    script:
    """
    iqtree -T 4 -s concat_all.txt -m MFP
    """
}
