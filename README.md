# Viral and Plasmid Sequence Recovery, Filtering, and Clustering Workflow

This repository documents the bioinformatic workflow used to recover **uncultivated viral genomes (UViGs)** and **plasmid sequences** from assembled metagenomes, perform quality filtering, cluster sequences into species-level **viral operational taxonomic units (vOTUs)** and **plasmid taxonomic units (PTUs)**, and refine the clustering by (i) reassigning fragmented representatives and (ii) promoting circular genomes as cluster representatives.

> Note: geNomad recovers both viral and plasmid sequences in a single end-to-end run. Downstream steps (length filter, CheckV, clustering, refinement) are performed independently for each sequence type.

---

## Overview

```
Metagenomes / Metatranscriptome
    │
    ▼
[geNomad end-to-end]   ── mines viruses AND plasmids in a single run
    │
    ├──────────────── viruses ────────────────┐
    │                                          │
    │                                          ▼
    │                              [Length ≥ 1 kb + CheckV]
    │                                          │
    │                                          ▼
    │              [vclust: prefilter → align → Leiden cluster]
    │                  (ANI 95%, qcov 85%, res 1.0)
    │                                          │
    │                                          ▼
    │                         [BLASTn all-vs-all + anicalc.py]
    │                                          │
    │                                          ▼
    │                       [refine_votus.py → final vOTUs]
    │
    └──────────────── plasmids ───────────────┐
                                               │
                                               ▼
                                    [Length ≥ 2 kb]
                                               │
                                               ▼
                       [vclust: prefilter → align → Leiden cluster]
                            (gANI 35%, res 0.9)
                                               │
                                               ▼
                               [BLASTn all-vs-all + anicalc.py]
                                               │
                                               ▼
                              [refine_ptus.py → final PTUs]
```

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| [geNomad](https://github.com/apcamargo/genomad) | v1.11.0 | Mining viral and plasmid sequences |
| geNomad database | v1.9 | Reference database |
| [CheckV](https://bitbucket.org/berkeleylab/checkv) | v1.0.1 | Viral genome quality assessment |
| CheckV database | v1.5 | Reference database |
| [vclust](https://github.com/refresh-bio/vclust) | v1.3.1 | ANI-based clustering |
| BLAST+ | ≥ 2.13 | All-vs-all nucleotide alignment |
| [anicalc.py](https://bitbucket.org/berkeleylab/checkv/src/master/scripts/anicalc.py) | from CheckV repo | ANI estimation |
| Python | ≥ 3.9 | Refinement scripts |
| pandas | ≥ 1.5 | Data handling |
| seqkit | ≥ 2.10 | Multifasta handling |

Install via the provided conda environment:

```bash
conda env create -f environment.yml
conda activate vp-workflow
```

---

## Step 1 — Mining viral and plasmid sequences (geNomad)

geNomad classifies and recovers **both viral and plasmid sequences** simultaneously in a single end-to-end run:

```bash
genomad end-to-end \
    --enable-score-calibration \
    --conservative \
    --lenient-taxonomy \
    --full-ictv-lineage \
    --sensitivity 7.0 \
    --cleanup \
    metagenome.fa \
    genomad_out \
    /path/to/genomad_db
```

**Parameter notes:**
- `--enable-score-calibration` — calibrates marker scores for more reliable classification.
- `--conservative` — uses stricter thresholds to minimize false positives.
- `--lenient-taxonomy` — allows partial taxonomic assignments down the lineage.
- `--full-ictv-lineage` — outputs full ICTV taxonomic ranks for viral hits.
- `--sensitivity 7.0` — increases MMseqs2 sensitivity for marker gene detection.
- `--cleanup` — removes intermediate files at the end of the run.

Output files used downstream:

- `genomad_out/<sample>_summary/<sample>_virus.fna` — viral sequences
- `genomad_out/<sample>_summary/<sample>_virus_summary.tsv` — viral metadata (topology, length, taxonomy, scores)
- `genomad_out/<sample>_summary/<sample>_plasmid.fna` — plasmid sequences
- `genomad_out/<sample>_summary/<sample>_plasmid_summary.tsv` — plasmid metadata (topology, length, scores)

After concatenating across samples, summary files are merged into:

- `all_virus_summaries.tsv`
- `all_plasmid_summaries.tsv`

---

## Step 2 — Quality filtering

### Viruses (length ≥ 1 kb + CheckV)

```bash
seqkit seq -m 1000 all_viruses.fna > viruses_1kb.fna

checkv end_to_end \
    viruses_1kb.fna \
    checkv_out \
    -t 16 \
    -d /path/to/checkv-db-v1.5
```

### Plasmids (length ≥ 2 kb)

```bash
seqkit seq -m 2000 all_plasmids.fna > plasmids_2kb.fna
```

---

## Step 3 — Species-level clustering with vclust

Following [MIUViG guidelines](https://www.nature.com/articles/nbt.4306), sequences are clustered using a three-step vclust pipeline (`prefilter` → `align` → `cluster`).

### vOTUs (ANI ≥ 95%, qcov ≥ 85%, Leiden resolution = 1)

```bash
# 1) Prefilter — keep only candidate pairs sharing ≥95% k-mer identity
vclust prefilter \
    -i viruses_1kb.fna \
    -o vOTUs/fltr.txt \
    --min-ident 0.95

# 2) Align — compute pairwise ANI between candidate pairs
vclust align \
    -i viruses_1kb.fna \
    -o vOTUs/ani.tsv \
    --filter vOTUs/fltr.txt

# 3) Cluster — Leiden clustering on ANI graph
vclust cluster \
    -i vOTUs/ani.tsv \
    -o vOTUs/clusters.tsv \
    --ids vOTUs/ani.ids.tsv \
    --algorithm leiden \
    --metric ani \
    --ani 0.95 \
    --qcov 0.85 \
    --leiden-resolution 1 \
    --out-repr
```

### PTUs (gANI ≥ 35%, Leiden resolution = 0.9)

```bash
# 1) Prefilter — looser thresholds reflecting greater plasmid sequence diversity
vclust prefilter \
    -i plasmids_2kb.fna \
    -o PTUs/fltr.txt \
    --min-kmers 20 \
    --min-ident 0.5

# 2) Align — compute pairwise ANI, applying ANI/qcov filters at output
vclust align \
    -i plasmids_2kb.fna \
    -o PTUs/ani.tsv \
    --filter PTUs/fltr.txt \
    --out-ani 0.70 \
    --out-qcov 0.50

# 3) Cluster — Leiden clustering using global ANI (gANI)
vclust cluster \
    -i PTUs/ani.tsv \
    -o PTUs/clusters.tsv \
    --ids PTUs/ani.ids.tsv \
    --algorithm leiden \
    --metric gani \
    --gani 0.35 \
    --leiden-resolution 0.9 \
    --out-repr
```

The `--out-repr` flag instructs vclust to designate one representative sequence per cluster.

---

## Step 4 — All-vs-all ANI validation

To detect cases where shorter representatives are actually fragments of longer genomes, an all-vs-all BLASTn search is performed on the cluster representatives, followed by ANI estimation with `anicalc.py`:

```bash
makeblastdb -in representatives.fna -dbtype nucl -out reps_db

blastn -query representatives.fna -db reps_db \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -out blast.tsv \
    -num_threads 16

python anicalc.py -i blast.tsv -o allVSall_ani.tsv
```

Run this independently for vOTU and PTU representatives.

---

## Step 5 — Refining cluster representatives

Two custom Python scripts (`refine_votus.py`, `refine_ptus.py`) apply two refinement rules:

1. **Fragment reassignment** — A representative `q` is reassigned to a longer representative `t` if:
   - `pid ≥ 95.0`
   - `qcov ≥ 85.0`
   - `tcov < 70` (i.e., `q` is largely contained within `t`, but `t` is much longer)
   - Transitive chains (A → B → C) are resolved to terminal targets.
2. **Circular trumps linear** — If a cluster representative is linear but contains circular members (geNomad topology = `DTR` or `ITR`), the longest circular member is promoted to representative.

Run as:

```bash
python workflow/05_refinement/refine_votus.py \
    --clusters vOTUs/clusters.tsv \
    --ani      vOTUs/allVSall_ani.tsv \
    --genomad  vOTUs/all_virus_summaries.tsv \
    --out      vOTUs/votu_clusters_updated_final.tsv

python workflow/05_refinement/refine_ptus.py \
    --clusters PTUs/clusters.tsv \
    --ani      PTUs/allVSall_ani.tsv \
    --genomad  PTUs/all_plasmid_summaries.tsv \
    --out      PTUs/ptu_clusters_updated_final.tsv
```

---

## Citation

If you use this workflow, please cite: TBD
