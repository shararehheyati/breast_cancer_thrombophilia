# Instructions for Downloading Raw Data

Due to GitHub's file size limits, the raw expression data cannot be stored directly in this repository.

## Steps to get the data:

1. Go to the [UCSC Xena Browser TCGA BRCA Dataset](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN&host=https%3A%2F%2Ftcga.xenahubs.net)

2. Click the "Download" button to get the file `TCGA.BRCA.sampleMap_HiSeqV2_PANCAN`

3. Rename the downloaded file to: `TCGA_BRCA_RNAseq_Expression.txt`

4. Place it in the `data/raw/` directory of your local repository clone.

## File Details:
- **File name:** `TCGA_BRCA_RNAseq_Expression.txt`
- **Size:** ~200-300 MB
- **Format:** Tab-separated values (.txt)
- **Content:** Gene expression matrix for TCGA Breast Cancer samples
