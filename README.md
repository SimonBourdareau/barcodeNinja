# barcodeNinja

Flexible demultiplexing of combinatorial barcodes in FASTQ files.


Simon Bourdareau, PhD
Zeitlinger Lab, Stowers Institute for Medical Research.

---

## What it is for

Standard demultiplexers assume barcodes live in the index reads, or at the very start of read 1. barcodeNinja assumes nothing: you describe where the barcodes are with a small design language, and it extracts, matches, counts, trims and splits accordingly.

It is the right tool when your barcodes are **at known offsets** inside reads or index files, and the structure is **combinatorial** enough that `bcl-convert` or a fixed-layout demultiplexer cannot express it.

A single read can carry a PCR handle you want ignored, a fixed barcode with several accepted variants, a UMI, and a lookup barcode — and barcodeNinja handles all four in one pass:

```
14X|GATC*ATGC*TTTT|8I|10L
```

It is probably **not** the right tool if you need position-aware deduplication (use Picard or umi_tools after alignment), cell calling and UMI error correction (use Cell Ranger), or long-read barcode detection (barcodes must be at fixed offsets).

---

## Building

Two backends. Both produce the same FASTQ output; they differ only in whether statistics can be written to an RDS file.

### Option A — pure C++

If you do not need the optional statistics in RDS format, no dependency beyond zlib is required:

```bash
make
```

Statistics are still available as a tab-separated file with `--statsTSV`, which works in every backend.

### Option B — with the libR library

If you want the RDS output, R must be available at compile time. The simplest route is a separate conda environment:

```bash
conda create -n R
conda activate R
conda install -c conda-forge r-base

# confirm R is discoverable
R CMD config --cppflags
R CMD config --ldflags
```

Then build from the directory containing `barcodeNinja.cpp`:

```bash
make clean
make BACKEND=R
conda deactivate
./barcodeNinja --version
```

`make BACKEND=R` bundles libR and the R runtime into `./lib`, so the resulting binary runs on nodes that have no R installed — which is why `conda deactivate` before running it is safe.

### Notes

Always `make clean` when switching backends. Object files live in `build/C/` and `build/R/`, so the two can coexist, but a stale binary at the top level will not be rebuilt otherwise.

`--version` reports which backend you built:

```
Compilation path :
libR, embedded R
R version at compile time: 4.5.3
```

or

```
Compilation path :
pure C++, no embedded R
```

---

## Quick start

Simplest useful case — an 8 bp sample barcode in the index file, matched against a lookup table, output split per sample:

```bash
barcodeNinja \
  -i run_I1.fastq.gz -1 "8L" \
  -r run_R1.fastq.gz \
  -l barcodes.tsv \
  -d "bcIndex1-t1" \
  -o results -p expt1
```

This produces `results/expt1_<BARCODE>_R1_demultiplexed.fastq.gz`, one file per barcode found, plus reads whose barcode matched nothing under the name `Un`.

---

## The design language

A **design** describes one read or index file as a sequence of tokens separated by `|`. Tokens are consumed left to right from the 5' end of the sequence.

There are four token types.

### `<n>L` — lookup

`8L` takes 8 bases and matches them against the lookup table supplied with `-l`. The assigned name goes into the read header.

Requires `-l`. barcodeNinja refuses to start if a design contains `L` and no lookup file is given.

### `<n>I` — include

`8I` takes 8 bases and accepts any sequence. Used for UMIs and random barcodes.

By default each distinct sequence is assigned a numeric ID (`ID1`, `ID2`, …). With `-u` the literal sequence is reported instead, which is what you normally want for a UMI you intend to use downstream.

### `<n>X` — exclude

`14X` takes 14 bases and discards them. The bases are consumed positionally — subsequent tokens start after them — but nothing is reported and no counter is created. Use this for PCR handles, spacers, and invariant sequence.

### Literal sequence — fixed barcode

`GATC` takes 4 bases and requires them to match. IUPAC ambiguity codes are supported: `GRTC` matches both `GATC` and `GGTC`.

Several accepted variants are separated by `*`:

```
GATC*ATGC*TTTT
```

All alternatives must be the same length. Reads that match none of them are flagged; with `-f` they are written to `*_rejected.fastq.gz` instead of the normal output.

### Putting it together

```
14X|GATC*ATGC*TTTT|8I|10L
```

reads as: skip 14, then a 4 bp fixed barcode with three accepted variants, then an 8 bp UMI, then a 10 bp lookup barcode. Total 36 bases consumed from the 5' end.

### Which option takes a design

| Option | Applies to |
|---|---|
| `-1, --bcIndex1` | index 1 file (`-i`), or the first index in the read header with `-H` |
| `-2, --bcIndex2` | index 2 file (`-I`), or the second header index with `-H` |
| `-3, --bcRead1` | read 1 (`-r`) |
| `-4, --bcRead2` | read 2 (`-R`) |

---

## Lookup table format

Tab-separated, two or three columns. Lines starting with `#` or `//` are ignored, as is anything after a `#` or `//` on a line.

```
# name        sequence     scope
BC01          ACGTACGT     SR|any
BC02          TTGGCCAA     SR|any
R1_A01        GGAATTCC     SR|bcIndex1-t1
R3_A01        CATGCATG     SR|bcIndex1-t5
```

**Column 1** — the name written into the read header and used for output file names.

**Column 2** — the barcode sequence. IUPAC codes allowed.

**Column 3** (optional, defaults to `SR|any`) — orientation and scope, separated by `|`.

*Orientation:*

- `S` — store the sequence as given
- `R` — store its reverse complement
- `SR` or `RS` — store both

*Scope:*

- `any` — this barcode may match at any `L` token
- `bcIndex1-t3` — only at token 3 of the index 1 design

Multiple scopes are separated by `*`: `SR|bcIndex1-t1*bcRead1-t2`.

Scoping matters for multi-round ligation barcoding, where the same sequence may be a valid Round 1 barcode but not a valid Round 3 barcode. **Token numbers count every token in the design, including `X` spacers.** In `8L|30X|8L|30X|8L` the three lookup barcodes are `t1`, `t3` and `t5`.

Identical names collapse into one counter, so several sequences can map to the same sample.

---

## Options

### Input

| Option | Meaning |
|---|---|
| `-r, --read1File` | read 1 FASTQ (gzipped). **Required.** |
| `-R, --read2File` | read 2, for paired-end |
| `-i, --index1File` | index 1 FASTQ |
| `-I, --index2File` | index 2 FASTQ |
| `-H, --bcIndexesInHeader` | read index sequences from the read headers instead of index files. Cannot be combined with `-i`/`-I`. |

### Barcodes

| Option | Default | Meaning |
|---|---|---|
| `-1 -2 -3 -4` | — | designs (see above) |
| `-l, --bclookupFilePath` | — | lookup table, required if any design uses `L` |
| `-m, --bcMaxMismatches` | 0 | mismatches allowed when matching a lookup barcode |
| `-t, --bcTrim` | off | remove all design tokens from the read sequence |
| `-u, --returnUMIasSequences` | off | report `I` tokens as sequence rather than numeric ID |
| `-f, --rejectNonMatchingFixedBarcodes` | off | send reads failing a fixed-barcode token to `*_rejected` |
| `-D, --deduplicateReads` | off | send exact duplicates to `*_duplicated` |

### Trimming

| Option | Default | Meaning |
|---|---|---|
| `-a, --adapterR1` | off | adapter to trim from read 1 |
| `-A, --adapterR2` | off | adapter to trim from read 2 |
| `--trimOverlapRequire` | 15 | minimum overlap to call an adapter |
| `--trimOverlapDiffLimit` | 2 | maximum mismatches in the overlap |
| `--trimdiffPercentLimit` | 0.1 | maximum mismatch fraction in the overlap |
| `--trimQuality` | 20 | 3' quality trimming threshold |
| `--trimMinStretchG` | 5 | minimum poly-G run to trim (two-colour chemistry) |
| `-k, --minLengthRead` | 22 | reads shorter than this after trimming are dropped. Minimum 12. |

### Output

| Option | Default | Meaning |
|---|---|---|
| `-o, --outputDir` | — | output directory. **Required.** Created recursively if absent. |
| `-p, --outputPrefix` | none | prefix for all output file names |
| `-d, --outputDemultiplexingOn` | — | tokens to split output files on |
| `--statsTSV` | off | write statistics to `<outputDir>/<prefix>_barcodeNinja_statistics.tsv`. Works in every backend. |
| `--checkDesign` | off | parse the designs, show how the first records would be split, then exit without writing anything |

### Other

| Option | Default | Meaning |
|---|---|---|
| `-c, --numThreads` | 1 | worker threads for gzip writing. 1 means sequential. |
| `--keepQuiet` | off | suppress the progress line |
| `--statsRDS` | off | write statistics to `<outputDir>/<prefix>_barcodeNinja_statistics.rds` (`BACKEND=R` only) |

---

## Output

### Read headers

The assigned barcode names are appended to the read name, `_` separated, in design order:

```
@rd0_BC01 1:N:0:0
```

With several tokens you get `@read_BC01_ACGTACGT_ID3`, in the order the tokens appear across index 1, index 2, read 1, read 2.

Reads rescued by mismatch tolerance carry an `r` prefix on the barcode name — `rBC01` rather than `BC01`. **If you parse cell barcodes out of headers, strip a leading `r`**, or the same barcode will appear as two.

Reads matching nothing are named `Un`.

### Files

Without `-d`:

```
<prefix>_R1_demultiplexed.fastq.gz
<prefix>_R2_demultiplexed.fastq.gz     (paired-end)
<prefix>_R1_rejected.fastq.gz          (with -f)
<prefix>_R1_duplicated.fastq.gz        (with -D)
```

With `-d "bcIndex1-t1"`, the barcode name is inserted:

```
<prefix>_BC01_R1_demultiplexed.fastq.gz
<prefix>_BC02_R1_demultiplexed.fastq.gz
...
<prefix>_Un_R1_demultiplexed.fastq.gz
```

Demultiplexing on several tokens produces one file per observed combination:
`-d "bcIndex1-t1|bcIndex2-t1"` gives `<prefix>_BC01_BC07_R1_demultiplexed.fastq.gz`.

### Statistics

Always printed to stdout at the end:

```
=== Read Statistics ===
Total reads processed: 4000000
Total reads written (passed length filter): 3987221
Unique reads: 3987221
```

**`--statsTSV`** writes the full tables to `<outputDir>/<prefix>_barcodeNinja_statistics.tsv` in every backend, including the plain `make` build. Long format, one row per record:

```
table	key	assigned_name	count	average_quality
Index_Sequence_1	CATGCATG	BC04	100	40,40,40,40,40,40,40,40
Read_Sequence_1	ACGTACGT	ID1	100	40,40,40,40,40,40,40,40
Combinations	BC04_ID1	NA	100	NA
Read1Length	48	NA	400	NA
ReadStatistics	TotalReads	NA	400	NA
```

In R:

```r
stats <- read.delim("results/expt1_barcodeNinja_statistics.tsv")
tables <- split(stats, stats$table)
tables$Combinations
```

**`--statsRDS`** (`BACKEND=R` only) writes the same information to `<outputDir>/<prefix>_barcodeNinja_statistics.rds` as a named list, readable with `readRDS()`:

| Element | Contents |
|---|---|
| `Index_Sequence_1`, `Index_Sequence_2`, … | per-token barcode tables: `Sequence`, `AssignedName`, `Count`, `AverageQuality` |
| `Read_Sequence_1`, … | same, for read-derived tokens |
| `Combinations` | every observed barcode combination and its count |
| `Read1Length`, `Read2Length` | read length distribution after trimming |
| `ReadStatistics` | totals |

`Count` columns are `numeric`, not `integer`, so counts above 2^31 are representable.

Both statistics files carry `--outputPrefix`, so several samples can be written to one output directory without overwriting each other.

## Worked examples

Read structures below are illustrative — **verify the offsets against your own library** before using them.

### ChIP-nexus, single-end

5 nt random barcode, 4 nt fixed sample barcode, insert, adapter read-through.

```bash
barcodeNinja \
  -r nexus_R1.fastq.gz \
  -3 "5I|CTGA*TGAC*GACT*ACTG" \
  -t -u -f \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -k 22 \
  -d "bcRead1-t2" \
  -o nexus_out -p expt1
```

`-t` strips all 9 nt so the output is pure insert. `-u` puts the random barcode into the header as sequence. `-f` sends reads failing the fixed barcode to `*_rejected`. Output splits four ways on the sample barcode.

Deduplication is deliberately left off: the canonical nexus dedup is UMI plus mapped 5' position, which can only be done after alignment.

### SHARE-seq

Three 8 nt ligation barcodes separated by 30 nt spacers.

```bash
barcodeNinja \
  -r share_R1.fastq.gz -R share_R2.fastq.gz \
  -i share_I1.fastq.gz \
  -1 "8L|30X|8L|30X|8L" \
  -3 "10I" \
  -l share_barcodes.tsv \
  -m 1 -u -t \
  -o share_out -p sample1
```

Scope each round's barcodes in the lookup table to `bcIndex1-t1`, `bcIndex1-t3`, `bcIndex1-t5` so a Round 1 sequence cannot match at Round 3.

**Do not add `-d` here.** 96³ is 884,736 combinations and therefore 884,736 simultaneously open gzip streams. The cell barcode goes into the header instead; split downstream.

### CUT&RUN / CUT&Tag, dual index

```bash
barcodeNinja \
  -i I1.fastq.gz -1 "8L" \
  -I I2.fastq.gz -2 "8L" \
  -r R1.fastq.gz -R R2.fastq.gz \
  -l illumina_indexes.tsv \
  -m 1 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -k 25 \
  -d "bcIndex1-t1|bcIndex2-t1" \
  -o cnr_out -p expt
```

### Single-molecule footprinting

Bisulfite-converted reads: demultiplex only, touch nothing else.

```bash
barcodeNinja \
  -i smf_I1.fastq.gz -1 "8L" \
  -r smf_R1.fastq.gz -R smf_R2.fastq.gz \
  -l smf_barcodes.tsv \
  -d "bcIndex1-t1" \
  -k 30 \
  -o smf_out -p run1
```

No `-a`: C→T conversion inflates mismatches against the adapter and makes trimming unreliable. No `-D`: dedup must be post-alignment.

### Amplicon panel with duplex UMIs

```bash
barcodeNinja \
  -r amp_R1.fastq.gz -R amp_R2.fastq.gz \
  -3 "12I|GTGACTGGAGTTCAGACGTGT" \
  -4 "12I" \
  -t -u -D -f \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -k 40 \
  -o amp_out -p panel
```

The fixed linker after the UMI acts as a structural check — with `-f`, reads missing it go to `*_rejected` rather than being force-fit. This is the case where `-D` is most defensible, since amplicons start at fixed positions.

---

## Checking that a run is correct

**Check the design before running anything.** `--checkDesign` parses the designs, prints the layout with byte offsets, shows how the first few records would be split, and exits without writing output:

```bash
barcodeNinja -r reads_R1.fastq.gz -3 "14X|GATC*ATGC*TTTT|8I|8L" -l barcodes.tsv -o /tmp --checkDesign
```

```
bcRead1   (34 bases consumed)
    t1  bases 1-14  (14 bp)  exclude
    t2  bases 15-18  (4 bp)  fixed  [GATC, ATGC, TTTT]
    t3  bases 19-26  (8 bp)  include (UMI)
    t4  bases 27-34  (8 bp)  lookup

  @cx0 1:N:0:0
    read1   NNNNNNNNNNNNNNGATCTTTCCTCAACGTACGT...
      t1 X  NNNNNNNNNNNNNN
      t2 S  GATC   <- matches
      t3 I  TTTCCTCA
      t4 L  ACGTACGT
```

A wrong offset shows up immediately as `<- NO MATCH` and a visible frameshift in the extracted tokens. Reads shorter than the design are reported explicitly. This turns a 40-minute mistake into a 2-second one.

**Then look at the `Un` fraction.** If most reads land in `Un`, the design offsets are wrong. This is the single most useful diagnostic the tool gives you.

**Check the `Combinations` table** in the RDS. For a plate design you should see roughly even depth across wells; for single-cell data, the usual knee.

**Check the read length distribution** after trimming. A large spike at short lengths suggests adapter trimming is firing more than it should.

**Run a subset first.** `head -c 100000000` on a gzipped FASTQ gives a truncated but valid stream for a fast sanity check.

---

## Troubleshooting

**`Too many demultiplexed output files open at once`**
One gzip stream is held open per output file for the whole run, and you have exceeded the descriptor limit. Either raise it (`ulimit -n 65535`) or demultiplex on fewer tokens and split downstream.

**`An index2 design was provided but the read headers contain only one index sequence`**
`-H` with `-2` set, on a single-index run. Drop `-2`.

**`The fastq files 'read1File' and 'read2File' are not sorted the same way`**
R1 and R2 are not in the same order, or are from different runs.

**`Empty token in bcDesign`**
A doubled or trailing `|` in a design string.

**`All '*'-separated alternatives must have the same length`**
Fixed-barcode variants differ in length. All alternatives in one token must match the token width.

**`--minLengthRead cannot be shorter than 12 bp`**
Reads below 12 bp are not useful to any aligner. Raise the value.

**Progress line scrolls instead of updating**
Output is redirected to a file rather than a terminal. Expected; the line is written once per million reads in that mode.

**Everything lands in `Un`**
The design does not match the actual read structure. Print a few reads (`zcat file.fastq.gz | head -8`) and count bases against your design.

---

## Behaviour worth knowing

**Mismatch rescue renames barcodes.** With `-m 1`, a rescued barcode becomes `rBC01`. This is intentional — it keeps rescued reads distinguishable — but it means rescued and exact reads land in different output files under `-d`, and header parsers must strip the leading `r`.

**Ambiguous matches are not rejected.** If a sequence is equidistant from two lookup barcodes at the allowed mismatch level, one is chosen and it is not currently defined which. Use `-m 0` where this matters.

**Deduplication is exact and pre-alignment.** `-D` compares the concatenated barcode and read sequences. A single sequencing error makes a duplicate look unique, and it does not know about mapping position. For fragment-level duplicates, dedup after alignment.

**Deduplication runs before adapter trimming.** Duplicates differing only in read-through length are not collapsed.

**Read headers are rewritten** as `@name_BARCODES 1:N:0:0`. The original index and filter fields are not preserved. Downstream tools that parse Illumina headers strictly may object.

**`X` tokens are consumed but never reported** — they still count for token numbering in `-d` and in lookup scopes.

**Memory** scales with the number of distinct barcodes and, with `-D`, with the number of unique reads. A large deduplicated run can use tens of gigabytes.

**Linux only, by design.** The build uses `/proc/self/exe`, `sys/sysinfo.h`, `getrusage` with Linux units, and POSIX headers throughout. There are no platform conditionals.

---

## Changes in this version

Several of these alter output relative to v1.8.1. Runs made with the older version may need repeating.

- **Adapter trimming no longer truncates reads that contain no adapter.** Previously any read could be cut at whatever offset had the fewest mismatches. Reads will be longer than v1.8.1 produced whenever `-a`/`-A` was used.
- **IUPAC ambiguity codes now match.** Previously a code such as `R` in a lookup entry or design matched nothing, silently sending all affected reads to `Un` or to `*_rejected`.
- **Quality scores above Q50 are handled.** Previously they were read as −1, which truncated reads during quality trimming and corrupted `AverageQuality`.
- **Counters are 64-bit.** Read counts above 2.1 billion and quality sums above ~54 million observations per barcode no longer overflow. RDS `Count` columns are now `numeric` rather than `integer`.
- **Reads shorter than the design are rejected** rather than read past the end of the buffer.
- **Paired-end ordering is checked** even when no barcode design is given.
- **Error messages are reported.** Previously almost every error printed as "Unknown error occurred."
- **`--minLengthRead` has a floor of 12 bp**; `-k 0` is no longer accepted.
- Build now uses `-O2`.

New in this version:

- **`--statsTSV`** writes the full statistics tables to a tab-separated file in every backend. Previously the pure C++ build computed the per-barcode counts and discarded them, so only an R-linked build could see the `Un` fraction or the combination distribution.
- **`--checkDesign`** validates a design against the first records of a file and exits, without writing output.
- **Both statistics files now carry `--outputPrefix`.** Previously the RDS was always written as `barcodeNinja_statistics.rds`, so several samples written to one output directory would silently overwrite each other's statistics. Scripts that hardcode the old name must be updated when a prefix is used.

If you have archived RDS files, a quick check for negative `AverageQuality` values will identify runs affected by the two quality bugs.
