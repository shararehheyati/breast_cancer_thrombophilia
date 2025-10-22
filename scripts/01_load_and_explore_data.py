# scripts/01_load_and_explore_data.py
"""
First analysis script: Load and explore thrombophilia gene expression in breast cancer
"""

import pandas as pd
import os

print("🔬 Starting analysis: Thrombophilia Genes in Breast Cancer")
print("=" * 60)

# 1. Setup paths
print("\n📁 Setting up file paths...")
data_dir = "data"
raw_data_dir = os.path.join(data_dir, "raw")
gene_list_path = os.path.join(data_dir, "thrombophilia_genes.txt")

# 2. Load thrombophilia gene list
print("📖 Loading thrombophilia gene list...")
try:
    with open(gene_list_path, 'r') as f:
        thrombophilia_genes = [line.strip() for line in f if line.strip()]
    print(f"✅ Found {len(thrombophilia_genes)} thrombophilia genes:")
    for i, gene in enumerate(thrombophilia_genes, 1):
        print(f"   {i:2d}. {gene}")
except FileNotFoundError:
    print("❌ ERROR: Gene list file not found!")
    print(f"   Make sure this file exists: {gene_list_path}")
    exit()

# 3. Check if raw data file exists
expression_file_path = os.path.join(raw_data_dir, "TCGA_BRCA_RNAseq_Expression.txt")
print(f"\n📊 Checking for expression data file...")
print(f"   Looking for: {expression_file_path}")

if not os.path.exists(expression_file_path):
    print("❌ ERROR: Expression data file not found!")
    print("   Please follow the instructions in data/raw/DOWNLOAD_INSTRUCTIONS.md")
    print("   to download the data and place it in data/raw/")
    exit()
else:
    print("✅ Expression data file found!")

# 4. Load a small portion of the data to explore
print("\n📈 Loading expression data (first look)...")
print("   This might take a moment for large files...")

try:
    # First, let's just look at the first few rows to understand the structure
    with open(expression_file_path, 'r') as f:
        first_lines = [f.readline() for _ in range(5)]
    
    print("\n📋 File structure preview:")
    for i, line in enumerate(first_lines):
        print(f"   Line {i}: {line[:100]}..." if len(line) > 100 else f"   Line {i}: {line.strip()}")
    
    # Now load the full data
    print("\n🔄 Loading full dataset...")
    expression_df = pd.read_csv(expression_file_path, sep='\t', index_col=0)
    
    print("✅ Data loaded successfully!")
    print(f"   Dataset shape: {expression_df.shape} (genes × samples)")
    print(f"   Number of genes: {expression_df.shape[0]}")
    print(f"   Number of samples: {expression_df.shape[1]}")
    
    # 5. Check which of our thrombophilia genes are in the dataset
    print(f"\n🎯 Finding our thrombophilia genes in the dataset...")
    genes_found = [gene for gene in thrombophilia_genes if gene in expression_df.index]
    genes_not_found = [gene for gene in thrombophilia_genes if gene not in expression_df.index]
    
    print(f"✅ Found {len(genes_found)} genes in the dataset:")
    for gene in genes_found:
        print(f"   ✓ {gene}")
    
    if genes_not_found:
        print(f"❌ {len(genes_not_found)} genes not found in dataset:")
        for gene in genes_not_found:
            print(f"   ✗ {gene}")

    # 6. Basic statistics for our genes of interest
    if genes_found:
        print(f"\n📊 Basic statistics for thrombophilia genes:")
        thrombophilia_data = expression_df.loc[genes_found]
        
        print("   Gene expression summary:")
        print(f"   Mean expression across all samples:")
        for gene in genes_found:
            mean_expr = thrombophilia_data.loc[gene].mean()
            print(f"      {gene}: {mean_expr:.2f}")
            
        # Save the subset for future analysis
        output_dir = os.path.join(data_dir, "processed")
        os.makedirs(output_dir, exist_ok=True)
        
        output_file = os.path.join(output_dir, "thrombophilia_expression_subset.csv")
        thrombophilia_data.to_csv(output_file)
        print(f"\n💾 Saved thrombophilia gene expression data to: {output_file}")

except Exception as e:
    print(f"❌ ERROR loading data: {e}")
    print("   This might be due to file format issues or memory limitations.")

print("\n" + "=" * 60)
print("✅ Script finished! Next steps:")
print("   1. Check the output above for any errors")
print("   2. If data loaded successfully, check the processed/ folder")
print("   3. We can now proceed to statistical analysis!")
