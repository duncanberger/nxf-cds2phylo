

def helpMessage() {
	log.info"""
Usage:
Mandatory arguments:
Optional arguments:
Skip metrics:
    """.stripIndent()
}

params.index = file(params.input)
params.outdir = "results"
params.phylo_method = "fasttree"
params.snpsites = true
params.partition = false
params.prefix = "out"

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
    if (params.snpsites == true & params.partition == false) {
    concatx(trim.out.collect())
    snpsites(concatx.out) }
    else {}
    if (params.phylo_method == "fasttree" ) {
    fasttree(snpsites.out) }
    if (params.phylo_method == "iqtree" & params.partition == true) {
    iqtree_partition(trim.out.collect()) }
    if (params.phylo_method == "iqtree" & params.partition == false) {
    iqtree_supertree(snpsites.out) }
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
    seqkit grep -j 4 -r -n -p ${genes} $baseDir/subset.test2.fasta > ${genes}.fa
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
    path("${params.prefix}.aln")

    script:
    """
    python $baseDir/scripts/concat_aln.py --input *.trimal --output ${params.prefix}.aln
    """
}

process snpsites {
    cpus = 1
    tag "X"
    publishDir "$params.outdir/", mode: 'copy'

    input:
    path("${params.prefix}.aln")

    output:
    path("${params.prefix}.snpsites")

    script:
    """
    snp-sites ${params.prefix}.aln -o ${params.prefix}.snpsites
    """
}

process fasttree {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("${params.prefix}.snpsites")

    output:
    path("${params.prefix}.fasttree")

    script:
    """
    FastTreeMP -gtr -nt -log logfile < ${params.prefix}.snpsites > ${params.prefix}.fasttree
    """    
}

process iqtree_supertree {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("${params.prefix}.snpsites")

    output:
    path("${params.prefix}.iqtree_supertree.treefile")

    script:
    """
    iqtree -T 4 -s ${params.prefix}.snpsites --prefix ${params.prefix}.iqtree_supertree -m MFP+ASC
    """
}

process iqtree_partition {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("*.trimal")

    output:
    path("${params.prefix}.iqtree_partition.treefile")

    script:
    """
    mkdir -p temp 
    mv *.trimal temp/
    iqtree -T 4 -p temp/ --prefix ${params.prefix}.iqtree_partition -m MFP
    """
}
