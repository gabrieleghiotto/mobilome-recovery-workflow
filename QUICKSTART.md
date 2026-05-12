# Quick Start Guide

## Installation (5 minutes)

### 1. Create conda environment

```bash
cd /Users/Gabriele/Documents/POSTDOC/GAnDi/workflow
conda env create -f environment.yml
conda activate vp-workflow
```

### 2. Download required databases

```bash
# geNomad database (requires internet, ~20GB)
mkdir -p databases
genomad download-db databases

# CheckV database (optional, ~10GB)
checkv download_db databases/checkv-db
```

### 3. Configure paths

Edit `config/workflow_config.yaml`:

```yaml
genomad_db: "/Users/Gabriele/Documents/POSTDOC/GAnDi/workflow/databases/genomad-db-v1.9"
checkv_db: "/Users/Gabriele/Documents/POSTDOC/GAnDi/workflow/databases/checkv-db-v1.5"
input_dir: "data/metagenomes"
```

## Usage (1-2 hours)

### Option A: Run entire pipeline

```bash
python workflow_orchestrate.py
```

### Option B: Run step-by-step

```bash
# Step 1: Mine viral and plasmid genomes
bash 01_genomad.sh

# Step 2: Quality filtering (MODULAR - recommended)
# - Uses: checkv contamination, completeness, complete_genomes, quality_summary
# - Best for: Large datasets, detailed analysis
bash 02_quality_filtering.sh

# OR Option B2: CheckV end-to-end (simpler, faster)
# bash 02_quality_filtering_endtoend.sh

# Step 3: Clustering
bash 03_vclust_clustering.sh

# Step 4: ANI validation
bash 04_all_vs_all_ani.sh

# Step 5: Refinement
bash 05_refine_representatives.sh
```

### Option C: Run from specific step onwards

```bash
python workflow_orchestrate.py --step 3 4 5
```

## CheckV Modes

**Default (02_quality_filtering.sh)**: Modular steps
- `checkv contamination` - Screen for contaminants
- `checkv completeness` - Assess genome completeness
- `checkv complete_genomes` - Extract complete genomes
- `checkv quality_summary` - Generate combined metrics

**Alternative (02_quality_filtering_endtoend.sh)**: Single command
- `checkv end_to_end` - All-in-one analysis

See [CHECKV_OPTIONS.md](CHECKV_OPTIONS.md) for comparison.

## Inputs

Place your metagenome FASTA files in:

```
data/metagenomes/
├── sample1.fna
├── sample2.fna
└── sample3.fna
```

Supported formats: `.fna`, `.fa`, `.fasta`

## Outputs

Final results in:

```
data/final/
├── votu_clusters_updated_final.tsv     # vOTU assignments
├── virus_representatives.fna            # vOTU representative sequences
├── ptu_clusters_updated_final.tsv       # PTU assignments (if plasmids found)
└── plasmid_representatives.fna          # PTU representative sequences (if plasmids found)
```

## Troubleshooting

### "Command not found: genomad"

```bash
conda activate vp-workflow
```

### Database not found errors

Set database paths in `config/workflow_config.yaml` or environment:

```bash
export GENOMAD_DB=/path/to/genomad-db
export CHECKV_DB=/path/to/checkv-db
```

### Out of memory

Reduce number of threads or process input in smaller batches:

```bash
export THREADS=8  # Use 8 cores instead of 16
bash 02_quality_filtering.sh
```

### Intermediate files taking too much space

Delete `data/clusters/vOTUs/*.blast*` and `data/clusters/PTUs/*.blast*` after step 4:

```bash
rm data/clusters/*/reps_db*
rm data/clusters/*/blast*
```

## Performance

Typical runtimes on a modern system with 16 cores:

- geNomad (per sample): 5-15 minutes
- Quality filtering: 2-5 minutes
- vclust clustering: 10-30 minutes (depends on sequence count)
- ANI validation: 5-15 minutes
- Refinement: <1 minute

**Total:** ~1-2 hours for small datasets, more for large metagenome collections

## Citation

If you use this workflow, please cite:

- **geNomad**: Camargo et al. (2023) https://www.biorxiv.org/content/10.1101/2023.03.05.531206v1
- **CheckV**: Nayfach et al. (2021) https://doi.org/10.1038/s41587-020-00774-7
- **vclust**: https://github.com/refresh-bio/vclust
- **MIUViG guidelines**: Roux et al. (2019) https://www.nature.com/articles/nbt.4306

## Next Steps

1. Explore cluster assignments: `pandas read_csv('data/final/votu_clusters_updated_final.tsv')`
2. Perform taxonomic classification on representatives
3. Generate comparative genomics analyses
4. Annotate genes in representative genomes

---

For more details, see `README.md`
