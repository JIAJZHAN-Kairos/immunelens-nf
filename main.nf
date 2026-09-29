nextflow.enable.dsl = 2

def quoteArg(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process PREPARE_IMMUNELENS {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'software_versions.txt'

    input:
    val sample_count

    output:
    path 'immunelens-library.tar.gz', emit: library
    path 'regions.tsv', emit: regions
    path 'software_versions.txt', emit: versions

    script:
    '''
    prepare_immunelens.sh
    '''

    stub:
    '''
    tar -czf immunelens-library.tar.gz -T /dev/null
    printf 'genome\tlocus\tchromosome\tstart\tend\n' > regions.tsv
    printf 'STUB RUN: no biological analysis performed\n' > software_versions.txt
    '''
}

process EXTRACT_COVERAGE {
    tag { sample }
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample), val(bam), path(local_inputs), val(corrections)
    path regions

    output:
    tuple val(sample), path("${sample}.coverage"), val(corrections)

    script:
    def source = local_inputs ? local_inputs.find { it.name.endsWith('.bam') }.toString() : bam
    """
    extract_coverage.py --sample ${quoteArg(sample)} --bam ${quoteArg(source)} --original-bam ${quoteArg(bam)} --regions ${quoteArg(regions)}
    """

    stub:
    """
    mkdir ${quoteArg(sample + '.coverage')}
    printf '{"genome":"hg38","stub":true}\n' > ${quoteArg(sample + '.coverage/sample.json')}
    """
}

process FIT_IMMUNELENS {
    tag { sample }
    publishDir "${params.outdir}/samples", mode: 'copy'

    input:
    tuple val(sample), path(coverage), val(corrections)
    path library

    output:
    tuple val(sample), path("${sample}.immunelens")

    script:
    """
    tar -xzf ${quoteArg(library)}
    export R_LIBS_USER="\$PWD/immunelens-library"
    cat > corrections.json <<'IMMUNELENS_CORRECTIONS'
    ${groovy.json.JsonOutput.toJson(corrections)}
    IMMUNELENS_CORRECTIONS
    run_immunelens.R ${quoteArg(sample)} ${quoteArg(coverage)} corrections.json
    """

    stub:
    """
    mkdir ${quoteArg(sample + '.immunelens')}
    printf 'sample\tlocus\tgenome\tstatus\tcell_fraction\tcorrection\tmessage\n' > ${quoteArg(sample + '.immunelens/estimates.tsv')}
    for locus in TCRA TCRB TCRG IGH; do
        printf '%s\t%s\thg38\tstub\tNA\tunadjusted\tSTUB RUN\n' ${quoteArg(sample)} "\$locus" >> ${quoteArg(sample + '.immunelens/estimates.tsv')}
    done
    """
}

process SUMMARISE_COHORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path manifest, stageAs: 'input_manifest.csv'
    path results, stageAs: 'collected/*'

    output:
    path 'cohort_estimates.tsv'
    path 'cohort_estimates_wide.tsv'
    path 'cohort_qc.json'
    path 'samplesheet.csv'

    script:
    """
    summarise_cohort.py ${quoteArg(manifest)} collected
    """
}

workflow {
    if (params.help) {
        log.info """
        ImmuneLENS Nextflow pipeline
        nextflow run JIAJZHAN-Kairos/immunelens-nf --input samplesheet.csv -profile docker
        Required: --input CSV with sample,bam columns (local or s3://).
        Optional: --outdir output directory. By default, results are placed beside the input manifest.
        Seqera: select an AWS Batch compute environment and the seqera profile.
        """.stripIndent()
    } else {
        if (!params.input) error 'Provide --input: a CSV manifest with sample,bam columns.'
        def manifest = file(params.input, checkIfExists: true)
        def samples = Channel.of(manifest).splitCsv(header: true, strip: true).collect().map { rows ->
            if (!rows) error 'The input manifest has no samples.'
            def names = [] as Set
            def bams = [] as Set
            rows.collect { row ->
                if (!row.containsKey('sample') || !row.containsKey('bam')) {
                    error 'The input manifest must contain sample,bam columns.'
                }
                def sample = row.sample?.trim()
                def bam = row.bam?.trim()
                if (!sample || !(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                    error "Invalid sample ID: ${sample}. Use letters, digits, underscore, hyphen or dot."
                }
                if (!bam || !bam.toLowerCase().endsWith('.bam') || bam.contains('\n') || bam.contains('\r')) {
                    error "Sample ${sample}: bam must be a BAM path or s3:// URI."
                }
                if (bam.startsWith('s3://') && !(bam ==~ /s3:\/\/[^\/]+\/.+/)) {
                    error "Sample ${sample}: invalid S3 BAM URI."
                }
                def local_inputs = []
                if (!bam.startsWith('s3://')) {
                    if (bam.contains('://')) error "Sample ${sample}: only local BAMs and s3:// BAMs are supported."
                    bam = file(bam, checkIfExists: !workflow.stubRun).toAbsolutePath().toString()
                    def candidates = [bam + '.bai', bam.replaceFirst(/\.bam$/, '.bai'),
                                      bam + '.csi', bam.replaceFirst(/\.bam$/, '.csi')]
                    def index = candidates.find { file(it).exists() }
                    if (!index) error "Sample ${sample}: a local BAI/CSI index is required."
                    local_inputs = [file(bam), file(index)]
                }
                if (!names.add(sample)) error "Duplicate sample ID: ${sample}."
                if (!bams.add(bam)) error "Duplicate BAM path for sample ${sample}."
                def corrections = [:]
                def cnColumns = ['TCRA_cn', 'TCRB_cn', 'TCRG_cn', 'IGH_cn']
                if (row.containsKey('purity') || cnColumns.any { row.containsKey(it) }) {
                    def required = ['purity', 'genome'] + cnColumns
                    if (!required.every { row.containsKey(it) }) {
                        error "Sample ${sample}: correction requires purity,genome,TCRA_cn,TCRB_cn,TCRG_cn,IGH_cn columns."
                    }
                    if (!(row.genome in ['hg19', 'hg38'])) error "Sample ${sample}: genome must be hg19 or hg38."
                    corrections.genome = row.genome
                    required.findAll { it != 'genome' }.each { column ->
                        def text = row[column]?.trim()
                        if (column != 'purity' && text in ['', 'NA']) {
                            corrections[column] = null
                        } else {
                            def number
                            try { number = Double.parseDouble(text ?: '') }
                            catch (Exception ignored) { error "Sample ${sample}: invalid ${column}: ${text}." }
                            if (!Double.isFinite(number) || number < 0 || (column == 'purity' && number > 1)) {
                                error "Sample ${sample}: ${column} must be finite and nonnegative; purity must be within [0,1]."
                            }
                            corrections[column] = number
                        }
                    }
                    cnColumns.each { column -> corrections[column + '_status'] = row[column + '_status'] ?: 'user_supplied' }
                }
                tuple(sample, bam, local_inputs, corrections)
            }
        }
        PREPARE_IMMUNELENS(samples.map { it.size() })
        EXTRACT_COVERAGE(samples.flatMap { it }, PREPARE_IMMUNELENS.out.regions)
        FIT_IMMUNELENS(EXTRACT_COVERAGE.out, PREPARE_IMMUNELENS.out.library)
        SUMMARISE_COHORT(manifest, FIT_IMMUNELENS.out.map { sample, result -> result }.collect())
        log.info "ImmuneLENS output directory: ${params.outdir}"
    }
}
