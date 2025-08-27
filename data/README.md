# Data Access Information

Due to the large size of the datasets, all data files required for the TransCisPredict analysis pipeline and for applying the trained protein prediction models are hosted on Synapse. This directory contains information for accessing these datasets.

## Data Access via Synapse

**Synapse Project**: [SynID: 69052240](https://www.synapse.org/Synapse:syn69052240/wiki/)
**Download Instructions**: See sections below

## Dataset Contents

**Structure**:
```
TransCisPredict/
|-- bim/
|   |-- `all_chromosomes.bim`      # Complete variant set across all autosomes (chromosomes 1-22)
|   |-- `chr01.bim`                # Chromosome 1
|   |-- `chr02.bim`                # Chromosome 2
|   |-- `chr03.bim`                # Chromosome 3
|   |-- `chr04.bim`                # Chromosome 4
|   |-- `chr05.bim`                # Chromosome 5
|   |-- `chr06.bim`                # Chromosome 6
|   |-- `chr07.bim`                # Chromosome 7
|   |-- `chr08.bim`                # Chromosome 8
|   |-- `chr09.bim`                # Chromosome 9
|   |-- `chr10.bim`                # Chromosome 10
|   |-- `chr11.bim`                # Chromosome 11
|   |-- `chr12.bim`                # Chromosome 12
|   |-- `chr13.bim`                # Chromosome 13
|   |-- `chr14.bim`                # Chromosome 14
|   |-- `chr15.bim`                # Chromosome 15
|   |-- `chr16.bim`                # Chromosome 16
|   |-- `chr17.bim`                # Chromosome 17
|   |-- `chr18.bim`                # Chromosome 18
|   |-- `chr19.bim`                # Chromosome 19
|   |-- `chr20.bim`                # Chromosome 20
|   |-- `chr21.bim`                # Chromosome 21
|   |-- `chr22.bim`                # Chromosome 22
|-- weights/
|   |-- `A1BG.csv`
|   |-- `AAMDC.csv`
|   |-- `AARSD1.csv`
|   |-- ... (2,339 total files)
|-- `weights.tar.gz`            # Compressed archive of all weight files
|-- `README.md`
```

### Variant Annotation Files (`bim/`)
- Contains PLINK-format variant information files (`.bim`) providing genomic variant information
- **Genome build**: GRCh37/hg19 coordinates
- **Population**: UK Biobank White British ancestry participants

**File format** (standard PLINK `.bim` format):
```
CHR    SNP           CM    BP         A1    A2
1      rs3131962     0     756604     A     G
1      rs115991721   0     767096     G     A
1      rs12562034    0     768448     A     G
```

**Columns**:
- **CHR**: Chromosome number (1-22)
- **SNP**: SNP identifier (rsID format)
- **CM**: Genetic distance in centiMorgans (set to 0)
- **BP**: Base-pair position (GRCh37/hg19 coordinates)
- **A1**: Reference allele
- **A2**: Alternative allele

**Organization**:
- Individual chromosome files (chr01.bim - chr22.bim) for memory-efficient, chromosome-specific processing
- Complete genome-wide file (all_chromosomes.bim) for comprehensive analyses
- All variants sorted by genomic position within chromosomes

### Protein Prediction Weights (`weights/`)
- Contains posterior weights for protein expression prediction
- **Total files**: 2,339 `.csv` files (one per protein)
- **Naming convention**: `{PROTEIN}.csv`
- **Example**: `A1BG.csv`

### Bulk Download (`weights.tar.gz`)
- Compressed archive containing all 2,339 protein weight files
- Download this single file if you need all protein prediction weights

**File Format**: `.csv` files with two columns:
- `variant_id`: SNP identifier (rsID format)
- `{method}_weights`: Prediction weights estimated using one of four methods:
  - `bayes_r_weights` - BayesR method
  - `susie_weights` - SuSiE method
  - `lasso_weights` - Lasso method
  - `enet_weights` - Elastic-net method

**Example**:
```csv
"variant_id","bayes_r_weights"
"rs3131962",3.79313579658051e-06
"rs115991721",-0.00142259048919599
"rs12562034",-0.000133076431117467
"rs4040617",-6.22890042961895e-06
```

## Download and Usage Instructions

### Step 1: Data Access
1. Register for Synapse account at [https://www.synapse.org/](https://www.synapse.org/)
2. Request access to the TransCisPredict project: [SynID: 69052240](https://www.synapse.org/Synapse:syn69052240/wiki/)
3. Follow data use agreement requirements
4. Download required datasets to your local analysis environment

### Step 2: Data Usage
1. **Locate variants**: Use the `.bim` files to map between rsIDs and genomic positions for variant matching with your genotype data
2. **Apply weights**: Use prediction weights in conjunction with variant annotations for protein expression prediction
   - Download individual protein files from `weights/` folder for specific proteins
   - Download `weights.tar.gz` for all proteins (extract with: `tar -xzf weights.tar.gz`)

## Statistical Methods

The protein prediction weights were estimated using four different statistical approaches:

1. **BayesR**: Bayesian regression with mixture priors
2. **SuSiE**: Sum of Single Effects model for fine-mapping
3. **Lasso**: L1-regularized linear regression
4. **Elastic-net**: Combined L1 and L2 regularized regression

Each protein's weights file contains results from one of these methods, as indicated by the column name suffix.


## Citation and Data Use

When using these data files, please:
- Cite the associated TransCisPredict manuscript [FIXME]
- Acknowledge the UK Biobank resource
- Comply with UK Biobank data use agreements
- Restrict use to academic research purposes