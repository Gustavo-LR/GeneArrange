# GeneArrange

**Annotation-based comparison and visualization of gene organization in linear and circular DNA molecules.**

GeneArrange compares annotated **gene identity, order, orientation and extremity adjacencies** across DNA sequences, including organelle genomes, plasmids, bacterial chromosomes and individual nuclear chromosomes. It groups equivalent arrangements, quantifies changed gene junctions, searches for candidate single-block signed inversions, and generates publication-oriented SVG figures and tabular results.

GeneArrange compares **annotations**, not DNA alignments: it does **not** call base-pair-resolution structural variants, validate genome assemblies, detect unannotated rearrangements, or infer chromosome rearrangement histories. Use validated, homologous gene identifiers for biologically meaningful comparisons.

## Requirements

- Python 3.10 or newer
- [Biopython](https://biopython.org/)
- Two or more GenBank files (`.gb`, `.gbk`, `.gbf`, `.genbank`), with **one sequence record per file**
- Gene features with `/gene` annotations. Features with only `/locus_tag` are not recognized in this version.

Install the dependency:

```bash
python -m pip install biopython
```

Run `python genearrange.py --help` for all options. Input may be a single GenBank file, a directory searched recursively, or a ZIP archive; comparisons require at least two records/files.

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

`--start-gene` is **not permitted** in linear mode. Gene order is never rotated across physical chromosome ends. For consistent SVGs, GeneArrange chooses the lexicographically canonical direction between the annotated forward order and its whole-molecule reverse complement. This decision changes **display direction only**, not the inferred structure.

### Choose which annotated genes appear in the figures (`--plot`)

**The plotting filter is separate from the comparison filter.** All gene-order
comparisons, breakpoint distances, exact structural groups, similarity clusters,
and seriation remain based exclusively on eligible **single-copy, contiguous,
known-strand** marker genes, regardless of `--plot`.

| `--plot` value | Gene features shown in overview, collapsed, group, and compact SVGs |
| --- | --- |
| `default` (default) | Original behavior: retained unique-copy, non-spliced comparison markers, subject to the normal presence cutoff |
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
- At least **three eligible requested genes must be shared by every accession** for exact structural grouping. A subset that does not meet this requirement produces a clear error; select additional shared markers. Selected names not annotated anywhere are reported in `excluded_genes.tsv` / `gene_presence.tsv`.
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

GeneArrange creates `genearrange_results/` unless `-o` selects another directory. The main outputs are:

| File | Contents |
| --- | --- |
| `overview.svg`, `overview_compact_view.svg` | All accessions, gene-level and compact views |
| `overview_collapsed.svg`, `overview_collapsed_compact_view.svg` | One representative per visualization group |
| `individual_image_groups/*.svg` | Per-group figures, plus compact counterparts |
| `samples.tsv` | Input accession metadata, explicit topology and normalized display direction |
| `gene_order.tsv`, `retained_genes.tsv` | **Comparison markers only**, per accession and across the dataset |
| `plotted_gene_order.tsv` | Every **displayed occurrence**, including duplicate-copy index, feature coordinates, parts, strand, and whether its name is a comparison marker |
| `gene_presence.tsv`, `excluded_genes.tsv` | Marker presence and why candidates were removed |
| `reference_comparisons.tsv`, `pairwise_comparisons.tsv` | Changed signed adjacencies and conserved blocks |
| `changed_junctions.tsv`, `exact_inversions.tsv`, `pairwise_blocks.tsv` | Candidate rearrangement descriptions |
| `breakpoint_distance_matrix.tsv`, `breakpoint_distance_pairs.tsv` | All-to-all distance matrix and normalized distances |
| `cluster_membership.tsv`, `cluster_summary.tsv` | Visualization group membership, representatives and thresholds |
| `display_order.tsv`, `visualization_groups.tsv` | Reproducible SVG row order and grouping |
| `README_results.md` | Run-specific parameters, tables and interpretation caveats |

The diagrams are **schematic, not proportional to base-pair coordinates**. Standard views display named gene arrows; compact views merge adjacent same-strand features and show gene identity through colored/texture stripes. Each repeated copy occupies its actual relative place among neighboring annotations. Styles are consistent across all figures of one run. The palette supports 34 vector texture families combined with 12 colors, for 408 base combinations before density variations. Unknown-strand features use a neutral rectangle.


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
