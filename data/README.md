# Data Access Information

Due to the large size of the datasets (>20GB), all data files required for the TransCisPredict analysis pipeline and for applying the trained protein prediction models are hosted on Synapse. This directory contains information for accessing these datasets.

## Data Access via Synapse

**Synapse Project**: [PLACEHOLDER_SYNAPSE_LINK]  
**Access Requirements**: [Details to be added]  
**Download Instructions**: See sections below

## Dataset Contents

### Variant Annotation Files
**Directory**: `bim_files/`
- Contains PLINK-format variant information files (.bim) with 646,634 variants total
- **Genome build**: GRCh37/hg19 coordinates
- **Population**: UK Biobank White British ancestry participants

**Structure**:
```
bim_files/
├── all_chromosomes.bim          # Complete variant set
├── chr01.bim                    # Chromosome 1
├── chr02.bim                    # Chromosome 2
├── chr03.bim                    # Chromosome 3
├── chr04.bim                    # Chromosome 4
├── chr05.bim                    # Chromosome 5
├── chr06.bim                    # Chromosome 6
├── chr07.bim                    # Chromosome 7
├── chr08.bim                    # Chromosome 8
├── chr09.bim                    # Chromosome 9
├── chr10.bim                    # Chromosome 10
├── chr11.bim                    # Chromosome 11
├── chr12.bim                    # Chromosome 12
├── chr13.bim                    # Chromosome 13
├── chr14.bim                    # Chromosome 14
├── chr15.bim                    # Chromosome 15
├── chr16.bim                    # Chromosome 16
├── chr17.bim                    # Chromosome 17
├── chr18.bim                    # Chromosome 18
├── chr19.bim                    # Chromosome 19
├── chr20.bim                    # Chromosome 20
├── chr21.bim                    # Chromosome 21
├── chr22.bim                    # Chromosome 22
└── README.md                    # Documentation for variant files
```

**File format** (standard PLINK .bim with header):
```
CHR    SNP           CM    BP         A1    A2
1      rs3131962     0     756604     A     G
1      rs115991721   0     767096     G     A
1      rs12562034    0     768448     A     G
```

**Columns**:
1. Chromosome number (1-22)
2. SNP identifier (rsID format)
3. Genetic distance in centiMorgans (0)
4. Physical position (GRCh37/hg19)
5. Reference allele (A1)  
6. Alternative allele (A2)

**Organization**: 
- Individual chromosome files (chr01.bim - chr22.bim) for memory-efficient analysis
- Complete genome-wide file (all_chromosomes.bim) for comprehensive analysis
- All variants sorted by genomic position within chromosomes

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

Weight files (`weights/{protein_name}_posterior_weights.csv`) contain:
```csv
variant_id,weight
chr1:12345:A:G,0.0234
chr1:67890:T:C,-0.0156
```

- **variant_id**: Chromosome:position:reference_allele:alternative_allele format
- **weight**: Effect size for the alternative allele in normalized protein expression units (in the column name you can tell which method it is generated from)

## Applications

These data files enable:
- Protein expression prediction in external cohorts with genotype data
- Proteome-wide association studies (PWAS) without measured proteomic data
- Investigation of molecular mechanisms underlying complex traits
- Identification of potential therapeutic targets