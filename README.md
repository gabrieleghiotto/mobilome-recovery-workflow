# Viral and Plasmid Sequence Recovery, Filtering, and Clustering Workflow

This repository documents the bioinformatic workflow used to recover **uncultivated viral genomes (UViGs)** and **plasmid sequences** from assembled metagenomes, perform quality filtering, cluster sequences into species-level **viral operational taxonomic units (vOTUs)** and **plasmid taxonomic units (PTUs)**, and refine the clustering by (i) reassigning fragmented representatives and (ii) promoting circular genomes as cluster representatives.

> **Note:** geNomad recovers both viral and plasmid sequences in a single end-to-end run. Downstream steps (length filter, CheckV, clustering, refinement) are performed independently for each sequence type.

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

## Quick Start

```bash
# Create conda environment
conda env create -f environment.yml
conda activate vp-workflow

# Run step-by-step
bash 01_genomad.sh
bash 02_quality_filtering.sh
bash 03_vclust_clustering.sh
bash 04_all_vs_all_ani.sh
bash 05_refine_representatives.sh
```

Or run all at once with the orchestration script:

```bash
python workflow_orchestrate.py --help
```

---

## Step-by-Step Documentation

### Step 1 — Mining viral and plasmid sequences (geNomad)

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

**Output files used downstream:**
- `genomad_out/<sample>_summary/<sample>_virus.fna` — viral sequences
- `genomad_out/<sample>_summary/<sample>_virus_summary.tsv` — viral metadata
- `genomad_out/<sample>_summary/<sample>_plasmid.fna` — plasmid sequences
- `genomad_out/<sample>_summary/<sample>_plasmid_summary.tsv` — plasmid metadata

See: `01_genomad.sh`

---

### Step 2 — Quality filtering

#### Viruses (length ≥ 1 kb + CheckV)

```bash
seqkit seq -m 1000 all_viruses.fna > viruses_1kb.fna

checkv end_to_end \
    viruses_1kb.fna \
    checkv_out \
    -t 16 \
    -d /path/to/checkv-db-v1.5
```

#### Plasmids (length ≥ 2 kb)

```bash
seqkit seq -m 2000 all_plasmids.fna > plasmids_2kb.fna
```

See: `02_quality_filtering.sh`

---

### Step 3 — Species-level clustering with vclust

Following [MIUViG guidelines](https://www.nature.com/articles/nbt.4306), sequences are clustered using a three-step vclust pipeline.

#### vOTUs (ANI ≥ 95%, qcov ≥ 85%)

```bash
vclust prefilter -i viruses_1kb.fna -o vOTUs/fltr.txt --min-ident 0.95
vclust align -i viruses_1kb.fna -o vOTUs/ani.tsv --filter vOTUs/fltr.txt
vclust cluster -i vOTUs/ani.tsv -o vOTUs/clusters.tsv \
    --algorithm leiden --metric ani --ani 0.95 --qcov 0.85 \
    --leiden-resolution 1 --out-repr
```

#### PTUs (gANI ≥ 35%)

```bash
vclust prefilter -i plasmids_2kb.fna -o PTUs/fltr.txt \
    --min-kmers 20 --min-ident 0.5
vclust align -i plasmids_2kb.fna -o PTUs/ani.tsv \
    --filter PTUs/fltr.txt --out-ani 0.70 --out-qcov 0.50
vclust cluster -i PTUs/ani.tsv -o PTUs/clusters.tsv \
    --algorithm leiden --metric gani --gani 0.35 \
    --leiden-resolution 0.9 --out-repr
```

See: `03_vclust_clustering.sh`

---

### Step 4 — All-vs-all ANI validation

To detect fragmented representatives, an all-vs-all BLASTn search is performed:

```bash
makeblastdb -in representatives.fna -dbtype nucl -out reps_db
blastn -query representatives.fna -db reps_db \
    -outfmt '6 std qlen slen' -max_target_seqs 10000 \
    -out blast.tsv -num_threads 16
python scripts/anicalc.py -i blast.tsv -o allVSall_ani.tsv
```

See: `04_all_vs_all_ani.sh`

---

### Step 5 — Refining cluster representatives

Two refinement rules are applied:

1. **Fragment reassignment** — A representative `q` is reassigned to a longer representative `t` if:
   - `pid ≥ 95.0`
   - `qcov ≥ 85.0`
   - `tcov < 70`
   - Transitive chains are resolved.

2. **Circular trumps linear** — If a cluster representative is linear but contains circular members, the longest circular member is promoted.

Run as:

```bash
python scripts/refine_votus.py \
    --clusters vOTUs/clusters.tsv \
    --ani vOTUs/allVSall_ani.tsv \
    --genomad vOTUs/all_virus_summaries.tsv \
    --out vOTUs/votu_clusters_updated_final.tsv

python scripts/refine_ptus.py \
    --clusters PTUs/clusters.tsv \
    --ani PTUs/allVSall_ani.tsv \
    --genomad PTUs/all_plasmid_summaries.tsv \
    --out PTUs/ptu_clusters_updated_final.tsv
```

See: `scripts/refine_votus.py`, `scripts/refine_ptus.py`

---

## Directory Structure

```
workflow/
├── README.md                                  # This file
├── environment.yml                            # Conda environment
├── workflow_orchestrate.py                    # Main orchestration script
│
├── 01_genomad.sh                             # Step 1: geNomad mining
├── 02_quality_filtering.sh                   # Step 2: Length & CheckV filtering
├── 03_vclust_clustering.sh                   # Step 3: Clustering
├── 04_all_vs_all_ani.sh                      # Step 4: ANI validation
├── 05_refine_representatives.sh              # Step 5: Refinement
│
├── scripts/
│   ├── anicalc.py                            # ANI calculation (from CheckV)
│   ├── refine_votus.py                       # vOTU refinement
│   ├── refine_ptus.py                        # PTU refinement
│   ├── merge_summaries.py                    # Merge geNomad outputs
│   └── extract_representatives.py            # Extract cluster representatives
│
├── config/
│   └── workflow_config.yaml                  # Configuration parameters
│
├── data/
│   ├── metagenomes/                          # INPUT: Raw metagenome sequences
│   ├── genomad/                              # geNomad results
│   ├── checkv/                               # CheckV results
│   ├── vclust/                               # vclust results
│   │   ├── vOTUs/
│   │   └── PTUs/
│   └── final/                                # FINAL OUTPUT
│       ├── votu_clusters_updated_final.tsv
│       ├── ptu_clusters_updated_final.tsv
│       ├── virus_representatives.fna
│       └── plasmid_representatives.fna
│
└── logs/                                     # Execution logs
```

---

## Usage Examples

### Run entire workflow

```bash
python workflow_orchestrate.py --input-dir data/metagenomes --output-dir data/final
```

### Run specific steps

```bash
# Just geNomad
bash 01_genomad.sh

# Just clustering (after geNomad)
bash 03_vclust_clustering.sh

# Just refinement (after clustering)
bash 05_refine_representatives.sh
```

### Custom configuration

Edit `config/workflow_config.yaml`:

```yaml
genomad:
  sensitivity: 7.0
  conservative: true

vclust:
  votu_ani: 0.95
  votu_qcov: 0.85
  ptu_gani: 0.35

filtering:
  virus_min_length: 1000
  plasmid_min_length: 2000
  checkv_enabled: true
```

Then run with custom config:

```bash
python workflow_orchestrate.py \
    --config config/workflow_config.yaml \
    --input-dir data/metagenomes
```

---

## Output Files

### Final cluster assignments

- `data/final/votu_clusters_updated_final.tsv` — vOTU cluster table
  - Columns: `representative_id`, `cluster_id`, `members`, `size`, `topology`, etc.

- `data/final/ptu_clusters_updated_final.tsv` — PTU cluster table
  - Columns: `representative_id`, `cluster_id`, `members`, `size`, etc.

### Representative sequences

- `data/final/virus_representatives.fna` — One representative per vOTU
- `data/final/plasmid_representatives.fna` — One representative per PTU

### Metadata

- `data/final/virus_summaries_annotated.tsv` — Viral metadata (length, topology, taxonomy, CheckV results)
- `data/final/plasmid_summaries_annotated.tsv` — Plasmid metadata

### Logs

- `logs/workflow.log` — Complete execution log with timing

---

## Citation

If you use this workflow, please cite: TBD

---

## References

- apcamargo/genomad: https://github.com/apcamargo/genomad
- CheckV: https://bitbucket.org/berkeleylab/checkv
- vclust: https://github.com/refresh-bio/vclust
- MIUViG: https://www.nature.com/articles/nbt.4306

---

**Last updated:** May 2026
**Contact:** Gabriele Ghiotto gabrieleghiotto@lbl.gov
