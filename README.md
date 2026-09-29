# ImmuneLENS Nextflow

Run [McGranahanLab/ImmuneLENS](https://github.com/McGranahanLab/ImmuneLENS)
on indexed WGS BAMs with **one required input: a CSV manifest**. Designed for
Seqera Platform with AWS Batch; local execution with Docker is also supported.

## Seqera launch

1. Add pipeline: `https://github.com/JIAJZHAN-Kairos/immunelens-nf`.
2. Select revision `v1.1.0` and profile `seqera`.
3. Select an existing Linux x86-64 AWS Batch compute environment with access to the BAM
   bucket and the manifest/output bucket. Its configuration supplies the queue,
   executor, work directory and AWS job role.
4. Enter just one pipeline parameter:

```yaml
input: s3://YOUR_BUCKET/immunelens/paad/samplesheet.csv
```

The Nextflow parameter schema generates the input form in Seqera.
If `outdir` is omitted, results go to
`s3://YOUR_BUCKET/immunelens/paad/results/immunelens-<UTC timestamp>/`.
Use `outdir` to choose a different destination.

Upload your existing manifest without changing its contents:

```bash
aws s3 cp /path/to/samplesheet.csv \
  s3://YOUR_BUCKET/immunelens/paad/samplesheet.csv --profile YOUR_PROFILE
```

The compute environment must permit `s3:GetObject` on BAMs, their BAI/CSI indexes
and the manifest, and the usual Nextflow read/write permissions on work and
output locations. Worker containers also need HTTPS access to GitHub, CRAN,
S3 and the public container registry. No AWS credentials are stored in this
repository. The default AWS region is `ap-southeast-2`.

**Fusion is not required.** The `seqera` profile explicitly disables Fusion
and runtime Wave. AWS Batch runs the fixed Docker image and Nextflow uses
ordinary S3 transfers for task inputs and outputs. BAMs are queried over HTTPS
after downloading their indexes, without a Fusion mount. The image hostname
`community.wave.seqera.io` is a public container registry; pulling this frozen
image does not require enabling Wave or Fusion in Seqera.

## Manifest

```csv
sample,bam
SAMPLE_A,s3://my-bucket/data/SAMPLE_A_tumor.bam
SAMPLE_B,s3://my-bucket/data/SAMPLE_B_tumor.bam
```

`sample` and `bam` are required and must be unique. IDs may contain letters,
digits, underscores, hyphens and dots, starting with a letter or digit.
Metadata columns other than the correction fields below are ignored.
All manifest rows are retained in the original
order. BAMs must be indexed and suitable for coordinate-based queries.
Indexes are discovered automatically as `file.bam.bai`, `file.bai`,
`file.bam.csi` or `file.csi` beside each BAM.

To apply purity/local-CN adjustment, include these additional columns in the
same manifest:

```csv
sample,bam,genome,purity,TCRA_cn,TCRB_cn,TCRG_cn,IGH_cn
SAMPLE_A,s3://my-bucket/data/SAMPLE_A_tumor.bam,hg38,0.5,3,3,3,3
SAMPLE_B,s3://my-bucket/data/SAMPLE_B_tumor.bam,hg38,0.6,2,2,NA,2
```

`purity` is a fraction within [0,1]. Each `*_cn` is absolute tumour copy
number for that locus, not global ploidy or a log ratio. CN values can be
noninteger and must be nonnegative. `genome` must match the detected BAM
build. `NA` or an empty CN skips adjustment for that locus and retains its
raw estimate. Optional `*_cn_status` columns explain missing or reviewed CNs.
The original two-column manifest remains supported for unadjusted analysis.

The publication uses nearby-gene CN proxies: TCRA—OR10G3, TCRB—PRSS58,
TCRG—STARD3NL and IGH—TMEM121. The preparation helper uses the same proxies
with PURPLE canonical gene CN calls and preserves the original sample order:

```bash
python3 bin/prepare_manifest.py --input samplesheet.csv \
  --purity-dir /path/to/purple_purity --gene-cn-dir /path/to/purple_cnv_gene \
  --output samplesheet.immunelens.csv
```

It creates a separate manifest, a per-locus source/checksum audit and a QC
JSON. `NO_TUMOR` inputs are rejected. Where a proxy gene's minimum and maximum
CN disagree, its CN is recorded as `NA` with `proxy_gene_heterogeneous` status;
no arbitrary average is used. The helper validates the hg38 gene annotation.

## Analysis

1. Validate the complete manifest before submitting sample jobs.
2. Install ImmuneLENS **1.0.3**, revision
   `df34702b46c9a8e0e201469e1e80e65c9996716d`, once per run. Download checksums
   are verified. Runtime R packages and samtools are recorded.
3. Identify hg19/hg38 and the `chr` prefix from chr7/chr14 lengths in each BAM
   header. Unknown builds fail explicitly. IGH is supported for hg38 only.
4. Read coverage only within upstream-defined TCRA, TCRB, TCRG and IGH
   analysis regions, including their normalization flanks. For S3 BAMs,
   download the small index and query the BAM over signed HTTPS range requests.
   Complete BAMs are not staged into Nextflow work directories. Signed URLs
   are kept out of scripts, output files and error logs.
5. Match upstream `getCovFromBam_WGS`: `samtools depth -q 20 -Q 20`, omitting
   uncovered positions. Use upstream model defaults, GC correction and flagged
   exon removal, including the default median coverage threshold of 15.
6. Export fractions, locus-specific segment usage, model fit, Shannon diversity
   and IGH class-switch metrics as provided by the upstream model.
   When correction columns are supplied, apply upstream `adjustImmuneLENS`
   to the summary, segment and model tables. Both raw and adjusted fractions
   are exported. Fractions exceeding [0,1] or the non-tumour fraction
   (`1-purity`) receive explicit QC status and are not clipped.
7. Gather exactly four locus records for every manifest sample. Missing or
   duplicate records fail the gather step. Low-coverage/no-estimate records are
   retained with `NA` fractions and explicit status; unexpected model failures
   fail the task rather than being silently converted into missing data.

The two-column manifest produces **unadjusted DNA-based estimates**.
The enriched manifest supplies purity and local CN for adjustment. Neither
format supplies matched-normal coverage for IGH germline/within-locus somatic
haplotype correction. Tumour-only IGH estimates and class-switch metrics
require particular caution at this polymorphic, copy-number-sensitive locus.
The three T-cell locus estimates are kept separate; no arbitrary consensus
fraction or immune hot/cold label is created. These estimates do not measure
immune function or spatial localization.

## Outputs

```text
<outdir>/
  cohort_estimates.tsv          # four rows per sample, including QC failures
  cohort_estimates_wide.tsv     # one row per sample, separate locus estimates/statuses
  cohort_qc.json                # manifest coverage, genome/status counts, correction state
  samplesheet.csv               # exact input manifest copy
  coverage/<sample>.coverage/
    TCRA.txt.gz, TCRB.txt.gz, TCRG.txt.gz, IGH.txt.gz
    bam_header.sam, sample.json
  samples/<sample>.immunelens/
    estimates.tsv, upstream_summaries.tsv
    <locus>.rds, <locus>.segments.tsv, <locus>.model.tsv
    <locus>.raw.rds             # original model when adjustment is applied
    sample.json, corrections.json, sessionInfo.txt
  pipeline_info/
    software_versions.txt, report.html, timeline.html, trace.tsv
```

Model files are present only when the upstream model returns an estimate.
IGH coverage is absent for hg19 and its status is `unsupported_genome`.
`cell_fraction` contains the adjusted estimate when valid correction inputs
are present, otherwise the raw estimate. `raw_cell_fraction` and
`adjusted_cell_fraction` are kept separately; inspect each locus's
`correction`, `cn_status` and `high_cell_fraction_flag`. The wide table has
matching per-locus columns and the cohort QC JSON records correction counts.
An execution success can include biological QC failures: inspect
`cohort_qc.json` and locus statuses before downstream analysis.

## Local execution and checks

```bash
nextflow run . --input /absolute/path/samplesheet.csv -profile docker

# Wiring test with two synthetic manifest rows; produces no biological estimates.
nextflow run . -profile test -stub-run

# Unit tests for genome detection, indexed coverage, and complete cohort gathering.
python3 -m unittest discover -s tests -v
```

Local BAMs and indexes are staged by Nextflow and mounted into Docker.
S3 access uses the worker AWS CLI credential
chain. On AWS Batch this is the job role; local Docker use must also expose
the selected credentials to the container.

Dependencies are a frozen public Seqera Containers image. Its build request is
in `conf/container-request.json`; ImmuneLENS is installed at runtime and is
not redistributed in that image or this repository. `restriktor` 0.6-50 is
also installed with a verified SHA-256 checksum. Reported provenance includes
the actual dependency versions, not only the pipeline version.

## License and citation

The wrapper code is MIT licensed. ImmuneLENS retains its
[academic-use license](https://github.com/McGranahanLab/ImmuneLENS/blob/main/LICENSE),
which includes redistribution restrictions. Users must comply with the
upstream terms; this repository does not grant rights to ImmuneLENS.

Please cite Bentham et al., *ImmuneLENS characterizes systemic immune
dysregulation in aging and cancer*, Nature Genetics (2025),
[doi:10.1038/s41588-025-02086-5](https://doi.org/10.1038/s41588-025-02086-5).
