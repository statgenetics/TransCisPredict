# Data Access Information

Due to the large size of the datasets (>20GB), all data files required for the TransCisPredict analysis pipeline and for applying the trained protein prediction models are hosted on Synapse. This directory contains information for accessing these datasets.

## Data Access via Synapse

**Synapse Project**: [PLACEHOLDER_SYNAPSE_LINK]  
**Access Requirements**: [Details to be added]  
**Download Instructions**: See sections below

## Dataset Contents

### Variant Annotation Files
**Directory**: `bim_files_per_LD_block/`
- Contains PLINK-format variant information files (.bim) organized by chromosome and LD blocks
- **Genome build**: GRCh37/hg19 coordinates

**Structure**:
```
bim_files_per_LD_block/
├── 01/
│   ├── 01_10583_1892607_unrelated_EUR_42644_individuals.bim
│   ├── 01_1892607_3582736_unrelated_EUR_42644_individuals.bim
│   └── ...
├── 02/
│   ├── 02_10133_1781022_unrelated_EUR_42644_individuals.bim
│   └── ...
├── 03/
└── ... (continues for all chromosomes)
```

**File format** (standard PLINK .bim):
```
1    rs3131962    0    756604     A    G
1    rs115991721  0    767096     G    A
1    rs12562034   0    768448     A    G
```

**Columns**:
1. Chromosome number (1-22)
2. Reference SNP ID (rs number)
3. Genetic distance in centiMorgans (0)
4. Physical position (GRCh37/hg19)
5. Reference allele (A1)  
6. Alternative allele (A2)

**File naming**: Files contain "EUR" indicating European ancestry subset, which consists of White British participants in this analysis.

### Prediction Weights
**Directory**: `weights/`
- Contains final trained protein prediction weights from Step 6 (Whole Sample Analysis)
- Incorporates both cis and trans effects across the genome
- Optimal statistical method selected per protein (BayesR, SuSiE, LASSO, or Elastic Net)
- Weight files for proteins that meet performance thresholds
- Each weight file contains variant IDs and their corresponding effect sizes for protein prediction

## Download and Usage Instructions

### Step 1: Data Access
1. Register for Synapse account at [https://www.synapse.org/](https://www.synapse.org/)
2. Request access to the TransCisPredict project: [PLACEHOLDER_SYNAPSE_LINK]
3. Follow data use agreement requirements
4. Download required datasets to your local analysis environment

### Step 2: Data Usage
1. **Locate variants**: Use the .bim files to map between rsIDs and genomic positions for variant matching with your genotype data
2. **Apply weights**: Use prediction weights in conjunction with variant annotations for protein expression prediction

## Technical Details

- **Coverage**: Genome-wide variants organized by LD blocks for both cis and trans effects
- **Quality Control**: Standard PLINK variant filtering applied
- **Coordinate System**: GRCh37/hg19 human genome reference
- **Ancestry**: UK Biobank White British ancestry participants

## Weight File Format

Weight files (`weights/{protein_name}_final_weights.csv`) contain:
```csv
variant_id,weight
chr1:12345:A:G,0.0234
chr1:67890:T:C,-0.0156
```

- **variant_id**: Chromosome:position:reference_allele:alternative_allele format
- **weight**: Effect size for the alternative allele in normalized protein expression units
- **Method**: Generated using optimal statistical method selected via cross-validation
- **Application**: Multiply genotype dosages by weights and sum to predict protein levels

## Applications

These data files enable:
- Protein expression prediction in external cohorts with genotype data
- Proteome-wide association studies (PWAS) without measured proteomic data
- Investigation of molecular mechanisms underlying complex traits
- Identification of potential therapeutic targets