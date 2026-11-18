# scbamqc

Alignment QC for single-cell BAM files, and cell selection from that QC.

- `scbamqc bam` writes one row of alignment statistics for each cell barcode.
- `scbamqc.qc` calculates a quality score for each barcode and selects cells.
- `scbamqc.pl` plots the statistics and the score.

## Install

```bash
pip install scbamqc
```

Add extras for more functions:

| Extra | Adds | Use it for |
|---|---|---|
| `bam` | pysam | `scbamqc bam` and `scbamqc tags` (Linux, macOS) |
| `anndata` | anndata, h5py | `scbamqc.io.read_10x_h5` |
| `plot` | matplotlib | `scbamqc.pl` and `scbamqc bam --plots` |
| `all` | all of the above | |

For example: `pip install "scbamqc[bam,plot]"`.

## Command line

The BAM file must be sorted by coordinate and indexed (`samtools index`).

```bash
# Cell Ranger
scbamqc bam possorted_genome_bam.bam -j 16 --out-prefix S1

# STARsolo SmartSeq: one read group for each cell
scbamqc bam Aligned.sortedByCoord.out.bam --barcode-tag RG
```

The tool writes `S1_per_cell.csv.gz`, with one row for each barcode. To also
write histograms to `S1_plots.png`, add `--plots`.

For all options, run `scbamqc bam --help`.

### Other pipelines

Different pipelines write the barcode and the mismatch count in different BAM
tags. Find the tags before you run `bam`:

```bash
scbamqc tags sample.bam
scbamqc bam sample.bam --barcode-tag CB --nm-tag NM
```

`bam` shows a warning if most reads do not have the tags you set.

> [!WARNING]
> Without a mismatch tag, `mean_nm` is empty (`NaN`). Then no barcode gets a
> quality score.

### Extra metrics

Add `--extra-metrics` to get more columns. The default columns do not change.
The new columns come after them. The tool writes the result to
`PREFIX_per_cell_extra.csv.gz`, not to `PREFIX_per_cell.csv.gz`. Thus the
output of a default run and an extra run can be in the same folder.

```bash
# Cell Ranger: the BAM tags give the read regions
scbamqc bam possorted_genome_bam.bam --extra-metrics

# Other aligners: a GTF file gives the read regions
scbamqc bam Aligned.sortedByCoord.out.bam --extra-metrics --gtf genes.gtf

# Smart-seq2: no strand
scbamqc bam sample.bam --barcode-tag RG --extra-metrics --gtf genes.gtf --strand unstranded
```

An alignment is one BAM record. A read with more than one alignment gives more
than one record. The denominator of each fraction is in the table.

| Column | Value | Denominator |
|---|---|---|
| `region_source` | Source of the region columns: `tags`, `gtf` or `none`. | |
| `frac_exonic` | Fraction of alignments in exons. | Alignments with a region. |
| `frac_intronic` | Fraction of alignments in introns. | Alignments with a region. |
| `frac_intergenic` | Fraction of alignments outside of genes. | Alignments with a region. |
| `frac_antisense` | Fraction of alignments on the opposite strand of a gene only. | Alignments with a region. |
| `n_umi` | Number of different UMIs on each contig, added for all contigs. | |
| `mean_reads_per_umi` | Primary alignments with a UMI, divided by `n_umi`. | `n_umi` |
| `frac_singleton_umi` | Fraction of UMIs with one primary alignment. | `n_umi` |
| `n_umi_counted` | Alignments with bit 8 of the Cell Ranger `xf` tag. | |
| `polya_tso_source` | Source of the polyA and TSO columns: `tags` or `seq`. | |
| `frac_polya` | Fraction of alignments with a polyA tail. | Alignments examined. |
| `frac_tso` | Fraction of alignments with the template switch oligo (TSO). | Alignments examined. |
| `frac_mapq_unique` | Fraction of alignments with MAPQ of `--mapq-unique` (default 255) or more. | All alignments (`n_reads`). |
| `frac_nm0` | Fraction of alignments with 0 mismatches (`--nm-tag`). | Alignments with the mismatch tag. |
| `frac_nm_ge2` | Fraction of alignments with 2 or more mismatches. | Alignments with the mismatch tag. |

The tool finds the sources from the first 20,000 mapped reads of the selected
contigs (`--contigs`, `--contig-prefix`). The sources are the same for all
barcodes.

#### Read regions

- The four region fractions add up to 1.
- If the BAM has the Cell Ranger `RE` tag, the tags give the region. An
  alignment with an `AN` tag and no `GX` or `TX` tag is antisense. Alignments
  without an `RE` tag are not in the denominator.
- If the BAM has no `RE` tag, the GTF file of `--gtf` gives the region. Only
  the aligned bases count. Soft clips, deletions (`D`) and splice gaps (`N`) do
  not count.
  - Exonic: 50% or more of the aligned bases are in exons.
  - Intronic: not exonic, and some bases are in an intron.
  - Antisense: bases are only in genes on the opposite strand.
  - Intergenic: all other alignments.
  - A gene without a strand (`.`) is on both strands.
  - A contig that is not in the GTF gives intergenic alignments. The log shows
    how many contigs are in the GTF.
- Without `RE` tags and without `--gtf`, the region columns are empty (`NaN`).
- To select the source, use `--region-source`:

  | `--region-source` | Source |
  |---|---|
  | `auto` (default) | The `RE` tags. If the BAM has no `RE` tag, the GTF. |
  | `tags` | The `RE` tags. If the BAM has no `RE` tag, the tool stops with an error. |
  | `gtf` | The GTF, also if the BAM has `RE` tags. You must also set `--gtf`. |

  To compare the GTF regions with the Cell Ranger regions, run the tool two
  times on the same BAM: one time with `tags` and one time with `gtf`.
- The two sources can give different results where genes on opposite strands
  overlap. The `RE` tag does not use the strand, so an alignment in an intron
  of its own gene and in an exon of the other gene can be exonic. The GTF
  source uses only the genes on the strand of the read, so it gives intronic.

Set the library strand with `--strand`. It applies to the GTF regions and to
the polyA and TSO sequence search.

| `--strand` | Strand of the transcript |
|---|---|
| `forward` (default) | Same strand as the read, for example 10x 3' R2. Read 2 of a pair has the opposite strand. |
| `reverse` | Opposite strand of the read. |
| `unstranded` | Both strands, for example Smart-seq2. No alignment is antisense. |

#### UMIs

- The UMIs come from the `--umi-tag` tag (default `UB`). Only primary
  alignments count.
- The UMIs are counted on each contig. A UMI with reads on two contigs counts
  as two UMIs.
- `n_umi` is not the UMI count of the count matrix. It counts all UMIs on the
  contig, also UMIs that Cell Ranger does not count, and does not use genes.
- `n_umi_counted` is for Cell Ranger BAMs. Cell Ranger sets bit 8 of the `xf`
  tag on one read of each molecule that it counts. With all contigs, the value
  should match the gene expression UMI count of the barcode in the Cell Ranger
  matrix. This is not yet verified on a real Cell Ranger BAM.
- If no read has the UMI tag, `n_umi`, `mean_reads_per_umi` and
  `frac_singleton_umi` are empty (`NaN`). This occurs for Smart-seq2 data.
- If no read has the `xf` tag, `n_umi_counted` is empty (`NaN`).

#### PolyA and TSO

- If the BAM has the Cell Ranger `pa` or `ts` tags, an alignment with a value
  more than 0 counts. All alignments are examined.
- If the BAM has no `pa` and no `ts` tag, the soft-clipped bases give the
  result. Alignments without a sequence (`SEQ` is `*`) are not examined.
  - PolyA: 10 or more `A` in a row in the clip at the 3' end of the
    transcript.
  - TSO: 12 or more bases of `AAGCAGTGGTATCAACGCAGAGTACATGGG` in the clip at
    the 5' end of the transcript.
  - The search assumes that the polyA tail and the TSO are on the transcript
    strand that `--strand` gives. With `unstranded`, the tool examines both
    strands.

#### Memory

With `--extra-metrics`, memory use also depends on the number of primary
alignments with a UMI on one contig. Each worker keeps 8 bytes for each of
these alignments. The measured peak of the UMI count for one contig is:

| Alignments per UMI | Peak for 10 million alignments |
|---|---|
| 1 | 0.25 GB |
| 3 | 0.14 GB |
| 5 | 0.12 GB |

Each worker has its own memory, so multiply by the number of workers (`-j`).

`qc_score` does not use the new columns.

## Python

```python
import scbamqc as sbq
from scbamqc.io import read_10x_h5

adata = read_10x_h5("raw_feature_bc_matrix.h5")
adata = sbq.add_bam_qc(adata, "S1_per_cell.csv.gz")
adata = sbq.basic_qc(adata)
sbq.qc_score(adata, min_reads=150)
sbq.select_cells(adata)

sbq.pl.qc_score_scatter(adata, min_score=0, min_counts=100)
sbq.pl.qc_score_hist(adata, threshold=0.9)

cells = adata[adata.obs["selected"] & (adata.obs["qc_score"] > 0.9)].copy()
```

| Step | Result |
|---|---|
| `add_bam_qc` | Keeps the barcodes in the BAM table. Adds the table columns to `obs`. |
| `basic_qc` | Removes barcodes with fewer than 5 genes and genes in fewer than 5 barcodes. Adds `total_counts`. |
| `qc_score` | Adds `obs["qc_score"]`, from 0 to 1. |
| `select_cells` | Adds `obs["gmm_label"]` and `obs["selected"]`. |

### The quality score

1. A barcode gets a score only if `n_reads > min_reads` and `mean_nm < 1`.
   Other barcodes get `NaN`.
2. A kernel density is fitted to three values of these barcodes:
   `mean_clip_frac`, `std_clip_frac` and `mean_base_qual` divided by its
   median.
3. The score is `log10(density + 1)`, scaled from 0 to 1.

A high score shows that the barcode aligns like most other barcodes.

Set `min_reads` for each library. The paper used values from 50 to 3,000.

`base_qual_scale` sets the divisor of `mean_base_qual`. The default is the
median of the scored barcodes. The divisor does not change the order of the
barcodes. It can change the score values only when the densities are near 1.

To use the same settings as the TgRN paper:

```python
sbq.qc_score(adata, min_reads=150, base_qual_scale=sbq.PAPER_BASE_QUAL_SCALE)
sbq.select_cells(adata, tol=sbq.PAPER_GMM_TOL)
```

### Cell selection

`select_cells` fits a Gaussian mixture with two components to
`log10(total_counts)` and `qc_score`. A barcode is selected if:

- it is in the component with the higher counts,
- `total_counts > min_counts` (default 100),
- `qc_score > min_score` (default 0).

The mixture fit runs until it converges. The paper stopped the fit early
(`tol=PAPER_GMM_TOL`). On the paper sample, the converged fit selects 15,536
cells and the paper fit selects 16,120. All cells of the converged fit are in
the paper set.

## R

The R package in [`r/`](r/README.md) does the same steps for Seurat objects.
It reads the table that `scbamqc bam` writes.

```r
remotes::install_github("umbibio/scbamqc", subdir = "r")
```

## Development

```bash
uv sync
uv run pytest
```

The tests in `tests/test_paper.py` compare the results with the paper. They
need the data in `../../data`, or in the folder set by `SCBAMQC_TEST_DATA`.
Without the data, pytest skips them.

One test in `tests/test_extra.py` reads the ToxoDB-68 *T. gondii* GTF. It
looks for `ToxoDB-68_TgondiiME49.gtf` in the same data folder. To use another
path, set `SCBAMQC_TOXODB_GTF`. Without the file, pytest skips the test.

## License

BSD 3-Clause.
