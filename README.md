# GeneArrange

**Annotation-based comparison and visualization of gene organization in linear and circular DNA molecules.**

GeneArrange is the general-purpose successor to **MitoArrange**. By default, it compares annotated **gene identity, order, orientation and extremity adjacencies** across DNA sequences, including organelle genomes, plasmids, bacterial chromosomes and individual nuclear chromosomes. It groups equivalent arrangements, quantifies changed gene junctions, searches for candidate single-block signed inversions, and generates publication-oriented SVG figures and tabular results. With `--no_cluster`, it instead **only visualizes annotated genes**, without requiring any genes shared among input records.

GeneArrange compares **annotations**, not DNA alignments: it does **not** call base-pair-resolution structural variants, validate genome assemblies, detect unannotated rearrangements, or infer chromosome rearrangement histories. Use validated, homologous gene identifiers for biologically meaningful comparisons.

## Requirements

- Python 3.10 or newer
- [Biopython](https://biopython.org/)
- Two or more GenBank files for structural comparison, or **one or more** for `--no_cluster` (`.gb`, `.gbk`, `.gbf`, `.genbank`), with **one sequence record per file**
- Gene features with `/gene` annotations. Features with only `/locus_tag` are not recognized in this version.

Install the dependency:

```bash
python -m pip install biopython
```

Run `python genearrange.py --help` for all options. Input may be a single GenBank file, a directory searched recursively, or a ZIP archive. The `--topology linear|circular` argument is **always required**.

## Quick start

### Circular sequences (mitochondria, chloroplasts, circular plasmids)

**Topology is mandatory.** This is intentionally independent of the GenBank `LOCUS` topology field.

```bash
python genearrange.py -i circular_genomes/ -o results_circular/ --topology circular
```

The default circular start/orientation gene is `cox1`. You can select another gene, even if it is not in the retained comparison set:

```bash
python genearrange.py -i circular_genomes/ -o results_circular/ --topology circular --start-gene rrn16 --reference Sample_A
```

### Linear sequences (bacterial or eukaryotic chromosomes, linear plasmids)

```bash
python genearrange.py -i chromosomes/ -o results_linear/ --topology linear
```

`--start-gene` is not permitted in **linear comparison mode**. Gene order is never rotated across physical chromosome ends. For consistent clustered SVGs, GeneArrange chooses the lexicographically canonical direction between the annotated forward order and its whole-molecule reverse complement. This decision changes **display direction only**, not the inferred structure. In linear **`--no_cluster`** mode, `--start-gene` is permitted as an **orientation cue**, but cannot move that gene to the beginning without misrepresenting the physical molecule.

### Visualization without structural grouping (`--no_cluster`)

Use `--no_cluster` when you want **only the gene-arrangement illustration** and no breakpoint distances, structural groups, representative selection, or similarity sorting. This is useful for poorly annotated genomes, records with fewer than three genes, or molecules with **no overlapping gene names**. One input record is sufficient.

```bash
python genearrange.py -i genomes/ -o results_visual/ --topology circular --no_cluster --plot all
```

By default, GeneArrange searches the **plotted gene names shared by every accession** for a gene with exactly one copy and a known strand in every record. It prioritizes `cox1` when available, otherwise takes the alphabetically first suitable name. The chosen gene defines a shared display origin and orientation for **circular molecules**; for **linear molecules**, it defines an orientation but does not rotate the physical ends. An explicit `--start-gene` overrides automatic anchor selection. If an anchor is unavailable or ambiguous in an accession, that accession follows its original GenBank orientation and the program emits a warning. If no suitable shared gene exists, **every accession is drawn in its original GenBank annotation order**.

```bash
python genearrange.py -i chromosomes/ -o results_visual/ --topology linear --no_cluster --start-gene dnaA --plot multicopy
```

In `--no_cluster`, `--plot default` means **single-copy, contiguous, known-strand genes within each accession**; unlike comparative mode, these genes do *not* need to be shared with other accessions or meet `--min_gene_presence`. `--plot multicopy`, `--plot all`, `--gene` and `--gene_file` also work, including the existing explicit multicopy override. When `--no_cluster` is present, grouping options and breakpoint thresholds are ignored with a warning. Rows follow alphabetical accession order, and the only figure is `overview.svg`.

### Choose which annotated genes appear in the figures (`--plot`)

**The plotting filter is separate from the comparison filter.** When comparison is enabled, all gene-order comparisons, breakpoint distances, exact structural groups, similarity clusters, and seriation remain based exclusively on eligible **single-copy, contiguous, known-strand** marker genes, regardless of `--plot`. When `--no_cluster` is selected, **no comparison is performed**.

| `--plot` value | Gene features shown in compact SVGs |
| --- | --- |
| `default` (default) | Comparison mode: retained unique-copy, non-spliced comparison markers, subject to the normal presence cutoff. `--no_cluster`: locally single-copy, contiguous, known-strand annotations, without a cross-sample presence requirement. |
| `multicopy` | All annotated non-spliced genes, **including every copy** of repeated gene names and genes below the presence cutoff. Two-part circular-origin joins are treated as contiguous in circular mode. |
| `all` | All named GenBank `gene` features, including repeated, multipart/spliced, and unknown-strand annotations. |

For example, to show a chromosome's duplicated genes without using them to infer structural groups:

```bash
python genearrange.py -i chromosomes/ -o results_multicopy/ --topology linear --plot multicopy
```

To draw every annotated named gene, including multipart genes:

```bash
python genearrange.py -i chromosomes/ -o results_all/ --topology linear --plot all
```

**Repeated copies are plotted individually in their genomic order**, interleaved
with the other genes using their original GenBank coordinates. They share a
color/texture identity and do not collapse into a single arbitrary-position
arrow. The same principle applies to the compact figures, where duplicated genes
appear as distinct colored stripes in their respective runs. Figures are ordered
by gene position but **not proportional to base-pair distances**.

In `--plot all`, a multipart gene is shown once at its first position in the
chosen orientation; the figure does not reconstruct the locations of all its
exons or trans-spliced segments. Unknown-strand features are drawn as rectangles
without an arrowhead. Gene features without a `/gene` identifier cannot be shown.

### Analyze only a subset of genes

Supply whitespace-separated gene identifiers or a comma-separated list.
`--gene` is an alias for `--genes`, and `--gene_file` is an alias for
`--genes-file`:

```bash
python genearrange.py -i chromosomes/ -o results_subset/ --topology linear --genes dnaA recA gyrA gyrB rpoB
```

```bash
python genearrange.py -i circular_genomes/ -o results_subset/ --topology circular --genes cox1,nad1,nad4,atp6,rrn16
```

Alternatively, create a text file containing gene names separated by newlines, commas, semicolons or whitespace:

```bash
python genearrange.py -i chromosomes/ -o results_subset/ --topology linear --genes-file selected_genes.txt
```

**Subset behavior:**

- Only requested gene names (case-insensitive) are retained as candidates; all other names are removed from the analysis **and plots**, including with `--plot all`.
- Every eligible selected marker is displayed wherever it is annotated, even when missing from other accessions. Missing genes remain available in the union legend.
- **Explicitly selecting a multicopy gene always plots every one of its copies**, even with `--plot default` and even if a particular copy is multipart. GeneArrange emits a warning and continues; repeated names remain **excluded from comparison**.
- When a gene subset is supplied, `--min_gene_presence` is ignored for comparison-marker eligibility. Multicopy, multipart (except a two-part circular-origin join), and unknown-strand names still cannot become comparison markers. `--plot all` may show excluded single-copy multipart features as well.
- **Exact structure (S) groups use the selected genes shared by all accessions** so that annotations missing in one accession do not force that accession into a separate group. Pairwise distances, similarity grouping and similarity display ordering project each pair onto its own shared annotated marker set.
- In comparison mode, at least **three eligible requested genes must be shared by every accession** for exact structural grouping. `--no_cluster` has **no minimum number of shared genes**. Selected names not annotated anywhere are reported in `statistics/excluded_genes.tsv` / `statistics/gene_presence.tsv` in comparison mode.
- Missing genes are **not counted as gene-order breakpoints**; gene content differences must be interpreted separately from comparisons of their shared gene order.

Without a subset argument, the original *comparison* selection behavior is preserved: all annotated genes are candidates and `--min_gene_presence 1.0` requires presence in every accession. Use `--min_gene_presence 0.8` to admit markers present in at least 80% of records. Non-default `--plot` values can display genes below that threshold without promoting them to structural markers.

For instance, to request four single-copy markers and a duplicated gene while keeping default plotting behavior:

```bash
python genearrange.py -i chromosomes/ -o results_subset/ --topology linear --plot default --gene dnaA recA gyrA rpoB rrn16
```

If `rrn16` has multiple annotated copies, **all copies are shown at their true gene-order positions**, but only eligible single-copy genes contribute to the comparisons. At least three qualifying shared markers are still required for subset grouping.

### Similarity-based grouping

```bash
python genearrange.py -i genomes.zip -o results/ --topology circular --group-by similarity --max-breakpoints 3 --collapsed_sub_all
```

Or use a normalized cutoff (range 0 to 1):

```bash
python genearrange.py -i chromosomes/ -o results/ --topology linear --group-by similarity --max-breakpoint-fraction 0.15 --display_order similarity
```

`--group-by` accepts `auto`, `folder`, `structure`, or `similarity`. `auto` selects folder grouping when at least two input folders exist, otherwise strict structural grouping. `similarity` uses deterministic complete-linkage clustering. `--display_order` can be `alphabetical` or `similarity` (default).

## How the topology-aware comparison works

### Circular mode

The retained signed gene order is oriented to the requested `--start-gene`, rotated to that gene, and compared using undirected gene **head/tail** extremity adjacencies around the circle. Rotation and whole-molecule reverse complementation are ignored by the structural signatures. The circular-origin multipart exception retains a *two-part gene* precisely crossing the physical zero/length boundary.

### Linear mode: a virtual neutral `END` marker

Linear DNA is represented by **all internal gene-to-gene junctions plus two virtual terminal connections**. Both chromosome extremities connect to the same orientation-neutral `END` identity, with no biological strand and **no physical gene added**. For illustration, `A+ — B- — C+` is internally represented as `END — A+ — B- — C+ — END`. The virtual marker never appears as an SVG arrow, and it is not included in the retained gene list.

This preserves **linear chromosome ends**, so `A—B—C—D` is **not** incorrectly equivalent to `C—D—A—B`, while both orientations of one molecule produce the **same adjacency signature**. Comparisons stay annotation-based; the virtual terminal markers correspond to the first and last **retained markers**, not known molecular telomere base pairs.

In both modes, gene directions (`+`/`-`) come from the GenBank feature strand. A whole-molecule reverse complement reverses gene order *and* flips every gene strand; signed head/tail adjacency sets are unchanged. The difference between topology modes is the presence of two neutral end connections versus the circular last-to-first connection.

## Marker selection and exclusions

For every candidate name, GeneArrange normally excludes it **from structural comparisons across the entire dataset** if it is:

- present as more than one `gene` feature in any record (multicopy);
- split among multiple segments (multipart), except for the precisely defined two-part circular-origin case in circular mode;
- annotated with undefined gene strand;
- below the requested presence fraction in default mode, or absent from every record in subset mode.

A marker set should consist of **orthologous genes with consistent annotation names**. The current tool does not run orthology inference, does not disambiguate paralog families, and does not use `/locus_tag` as a substitute for `/gene`. When examining large bacterial or nuclear chromosomes, this constraint is especially important.

## Outputs

**Compact SVGs are now the only figure format.** The former `_compact_view.svg` suffix has been removed because the compact presentation is the standard. The previous full-size, gene-labeled SVGs are no longer generated. Every `.tsv` report is placed in **`statistics/`**, rather than the top-level results directory.

### Standard comparison run

```text
results/
├── overview.svg                       # All accessions, compact style
├── overview_collapsed.svg             # One representative per group
├── individual_image_groups/
│   ├── S01.svg                        # Compact-only per-group figures
│   └── ...
├── statistics/
│   ├── samples.tsv
│   ├── retained_genes.tsv
│   ├── gene_order.tsv
│   ├── plotted_gene_order.tsv
│   ├── gene_presence.tsv
│   ├── excluded_genes.tsv
│   ├── reference_comparisons.tsv
│   ├── pairwise_comparisons.tsv
│   ├── pairwise_blocks.tsv
│   ├── changed_junctions.tsv
│   ├── exact_inversions.tsv
│   ├── all_pairs.tsv
│   ├── breakpoint_distance_matrix.tsv
│   ├── breakpoint_distance_pairs.tsv
│   ├── display_order.tsv
│   ├── visualization_groups.tsv
│   ├── cluster_membership.tsv
│   └── cluster_summary.tsv
└── README_results.md
```

### Visualization-only run (`--no_cluster`)

```text
results/
├── overview.svg                       # All accessions, compact style
├── statistics/
│   ├── samples.tsv                    # Accession and anchor metadata
│   ├── plotted_gene_order.tsv         # Every displayed occurrence and its coordinates
│   ├── gene_presence.tsv              # Gene presence and total annotated copies
│   └── display_order.tsv              # Figure row positions
└── README_results.md
```

No comparison-specific statistics, collapsed figures or individual-group figures are created in `--no_cluster` mode. Reusing an existing output directory automatically removes obsolete **GeneArrange-generated** SVGs and stale TSV reports from previous runs, including earlier `_compact_view.svg` names. Other unrelated files are not deleted.

**Figure interpretation:** compact SVGs merge adjacent same-strand gene features into arrows. Their internal colored and patterned stripes show gene identities; repeated genes retain their original genomic positions. Diagrams are schematics of gene order, **not proportional to base pairs**. Unknown-strand features use rectangles and multipart genes are shown once at a representative annotated position.

### Interpreting distances

`changed_adjacencies` counts signed extremity adjacencies seen in the reference or first sample that are absent from the other sample after projection onto pairwise shared genes; it is **not** the number of individual evolutionary events. Because both molecules have the same projected gene set, the distance is symmetric. `normalized_distance` divides by the number of common markers in circular mode, or **common markers + 1** in linear mode (accounting for both ends).

Exact structure groups compare signed adjacency signatures, but in subset mode are inferred only from genes shared across all samples. Similarity clusters use a user-controlled threshold and complete linkage; row ordering by similarity is just a visualization heuristic, **not a phylogeny**.

## Important limitations

1. GenBank inputs must each contain a **single DNA record**, and all records in a run must be treated as the same user-selected topology; the tool does not yet compare a multi-chromosome genome as one unit.
2. Comparisons require compatible gene names and comparable markers. A shared name does not prove orthology.
3. Intergenic-only changes, rearrangements lacking annotated marker boundaries, different nucleotide sequences with the same gene order, and heteroplasmic isoforms are not detected.
4. Missing annotations can hide genuine differences. The subset filter supports incomplete presence, but comparisons may become less informative as marker overlap decreases.
5. Repeated genes, typical intron-containing multipart gene annotations, trans-spliced genes and unknown-strand genes are excluded **from structural comparison**, but may be displayed by `--plot` or the explicit multicopy subset override. They can still carry uncertain orthology.
6. The exact signed-inversion search is computationally expensive for very large marker sets. SVG diagrams with hundreds or thousands of genes will also be difficult to interpret without further block-level visualization.
7. `--topology` is intentionally **mandatory** and overrides the topology recorded in the GenBank file. Choose the actual molecular topology carefully; the software does not validate it from sequence evidence.
8. Gene-array arrows display only named genes permitted by the requested plotting filter; differences in intergenic distances are not scaled or analyzed. Multipart genes are represented by one positional symbol, not exon-accurate geometry.
