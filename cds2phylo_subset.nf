def helpMessage() {
        log.info"""
Usage:
        nextflow run cds2phylo_subset.nf --input gene.list --fasta cds.fasta

Mandatory arguments:
        --input                         Path to file containing list of geneIDs to include
        --fasta                         Path to file containing fasta sequences to process
        --subset                        File of sample IDs to include

Optional arguments:
        --prefix                        Output directory ["prefix"]
        --phylo_model                   Use Fasttree or IQTREE ["fasttree"]
        --iqtree_parameters             Other parameters to pass to IQTREE [""]
        --partition                     Do partitioned analysis instead of supermatrix [false]
        --iqtree_model_supermatrix      IQTREE supermatrix analysis model ["MFP+ASC"]
        --iqtree_model_partition        IQTREE partioned analysis model ["MFP"]
        --skip_tree                     Skips creating the phylogeny [false]
    """.stripIndent()
}

if (!params.fasta){
    log.error "Error: '--fasta' parameter is missing"
    helpMessage()
    exit 0
}

if (!params.input){
    log.error "Error: '--input' parameter is missing"
    helpMessage()
    exit 0
}

params.index = file(params.input)
params.outdir = "results"
params.phylo_method = "fasttree"
params.snpsites = true
params.partition = false
params.subset = "input.list"
params.iqtree_model_supermatrix = "MFP+ASC"
params.iqtree_model_partition = "MFP"
params.iqtree_parameters = ""
params.skip_tree = false

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
    if (params.phylo_method == "fasttree" & params.skip_tree == false) {
    fasttree(snpsites.out) }
    if (params.phylo_method == "iqtree" & params.partition == true & params.skip_tree == false) {
    iqtree_partition(trim.out.collect()) }
    if (params.phylo_method == "iqtree" & params.partition == false & params.skip_tree == false) {
    iqtree_supertree(snpsites.out) }
    if (params.skip_tree == true) {}
}

process split {
    cpus = 1
    tag "X"

    input:
    tuple val(genes)

    output:
    tuple val(genes), path("${genes}.fa")

    script:
    """
    seqkit grep -j 1 -r -n -p ${genes} <(grep -v '^=' $baseDir/${params.fasta}) | seqkit grep -j1 -r -n -f $baseDir/${params.subset} > ${genes}.fa
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
    path("${params.subset}.aln")

    script:
    """
    python $baseDir/scripts/concat_aln.py --input *.trimal --output ${params.subset}.aln
    """
}

process snpsites {
    cpus = 1
    tag "X"
    publishDir "$params.outdir/", mode: 'copy'

    input:
    path("${params.subset}.aln")

    output:
    path("${params.subset}.snpsites")

    script:
    """
    snp-sites ${params.subset}.aln -o ${params.subset}.snpsites
    """
}

process fasttree {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("${params.subset}.snpsites")

    output:
    path("${params.subset}.fasttree")

    script:
    """
    FastTreeMP -gtr -nt -log logfile < ${params.subset}.snpsites > ${params.subset}.fasttree
    """
}

process iqtree_supertree {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("${params.subset}.snpsites")

    output:
    path("${params.subset}.iqtree_supertree.treefile")

    script:
    """
    iqtree -T 4 -s ${params.subset}.snpsites --subset ${params.subset}.iqtree_supertree -m ${params.iqtree_model_supermatrix} ${params.iqtree_parameters}
    """
}

process iqtree_partition {
    cpus = 4
    tag "X"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path("*.trimal")

    output:
    path("${params.subset}.iqtree_partition.treefile")

    script:
    """
    mkdir -p temp
    mv *.trimal temp/
    iqtree -T 4 -p temp/ --subset ${params.subset}.iqtree_partition -m ${params.iqtree_model_partition} ${params.iqtree_parameters}
    """
}
