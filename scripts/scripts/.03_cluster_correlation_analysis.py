# scripts/03_cluster_correlation_analysis.py
"""
آنالیزهای پیشرفته خوشه‌بندی و همبستگی برای ژن‌های ترومبوفیلیا
هدف: شناسایی الگوهای بیان و ارتباطات بین ژن‌ها
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.stats import spearmanr
import os

print("🔬 Starting Cluster & Correlation Analysis...")
print("=" * 60)

# 📚 تنظیمات حرفه‌ای
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# ایجاد پوشه‌های خروجی
os.makedirs('results/figures', exist_ok=True)
os.makedirs('results/tables', exist_ok=True)

# 📁 بارگذاری داده‌ها
print("\n📖 Loading processed data...")
try:
    expression_data = pd.read_csv("data/processed/thrombophilia_expression_subset.csv", index_col=0)
    print(f"✅ Data loaded! {expression_data.shape[0]} genes, {expression_data.shape[1]} samples")
except FileNotFoundError:
    print("❌ Processed data not found! Please run 01_load_and_explore_data.py first")
    exit()

# 🔥 ۱. Clustered Heatmap با خوشه‌بندی پیشرفته
print("\n🔥 ۱. Creating Advanced Clustered Heatmap...")

plt.figure(figsize=(16, 12))

# ایجاد clustered heatmap
clustered_grid = sns.clustermap(expression_data,
                               cmap='RdBu_r',
                               center=0,
                               metric='correlation',
                               method='average',
                               figsize=(16, 12),
                               yticklabels=True,
                               xticklabels=False,
                               dendrogram_ratio=(0.1, 0.1),
                               cbar_kws={'label': 'Expression Level'})

plt.suptitle('Clustered Heatmap: Thrombophilia Genes Expression Patterns\n(Genes Clustered by Correlation)', 
             fontsize=16, fontweight='bold', y=0.95)

# ذخیره نمودار
plt.savefig('results/figures/Figure5_clustered_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('results/figures/Figure5_clustered_heatmap.pdf', bbox_inches='tight')
print("✅ Figure 5: Clustered Heatmap saved!")

# 📊 ۲. Correlation Heatmap دقیق
print("\n📊 ۲. Creating Detailed Correlation Heatmap...")

# محاسبه ماتریس همبستگی
correlation_matrix = expression_data.T.corr(method='spearman')

plt.figure(figsize=(12, 10))

# ایجاد heatmap همبستگی
mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))

sns.heatmap(correlation_matrix,
           mask=mask,
           cmap='coolwarm',
           center=0,
           annot=True,
           fmt='.2f',
           square=True,
           cbar_kws={'label': 'Spearman Correlation Coefficient'})

plt.title('Correlation Matrix: Thrombophilia Genes Co-expression Patterns', 
          fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()

plt.savefig('results/figures/Figure6_correlation_heatmap.png', dpi=300, bbox_inches='tight')
print("✅ Figure 6: Correlation Heatmap saved!")

# 📈 ۳. Network Visualization از همبستگی‌های قوی
print("\n📈 ۳. Creating Strong Correlations Plot...")

# فیلتر کردن همبستگی‌های قوی
correlation_threshold = 0.6
strong_correlations = correlation_matrix[abs(correlation_matrix) > correlation_threshold]
strong_correlations = strong_correlations.fillna(0)

plt.figure(figsize=(10, 8))

sns.heatmap(strong_correlations,
           cmap='coolwarm',
           center=0,
           annot=True,
           fmt='.2f',
           square=True,
           cbar_kws={'label': f'Correlation > {correlation_threshold}'})

plt.title(f'Strong Correlations Between Thrombophilia Genes\n(Threshold: |r| > {correlation_threshold})', 
          fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()

plt.savefig('results/figures/Figure7_strong_correlations.png', dpi=300, bbox_inches='tight')
print("✅ Figure 7: Strong Correlations Heatmap saved!")

# 📋 ۴. آنالیز آماری همبستگی‌ها
print("\n📋 ۴. Performing Statistical Analysis of Correlations...")

# پیدا کردن قوی‌ترین همبستگی‌ها
correlation_pairs = []
for i in range(len(correlation_matrix.columns)):
    for j in range(i+1, len(correlation_matrix.columns)):
        gene1 = correlation_matrix.columns[i]
        gene2 = correlation_matrix.columns[j]
        corr_value = correlation_matrix.iloc[i, j]
        correlation_pairs.append({
            'Gene1': gene1,
            'Gene2': gene2,
            'Correlation': corr_value,
            'Abs_Correlation': abs(corr_value)
        })

# تبدیل به DataFrame و مرتب‌سازی
correlation_df = pd.DataFrame(correlation_pairs)
correlation_df = correlation_df.sort_values('Abs_Correlation', ascending=False)

print("\n🔍 Top 10 Strongest Gene Correlations:")
print(correlation_df.head(10).round(3))

# 💾 ۵. ذخیره نتایج برای مقاله
print("\n💾 ۵. Saving Analysis Results for Paper...")

correlation_matrix.to_csv('results/tables/correlation_matrix.csv')
correlation_df.to_csv('results/tables/top_correlations.csv', index=False)

with pd.ExcelWriter('results/tables/cluster_correlation_results.xlsx') as writer:
    correlation_matrix.to_excel(writer, sheet_name='Correlation_Matrix')
    correlation_df.head(20).to_excel(writer, sheet_name='Top_Correlations')
    expression_data.describe().to_excel(writer, sheet_name='Expression_Stats')

print("✅ All results saved to 'results/tables/'")

print("\n" + "=" * 60)
print("🎉 CLUSTER & CORRELATION ANALYSIS COMPLETED!")
print("📊 Generated 3 advanced figures:")
print("   🔥 Figure 5: Clustered Heatmap with dendrograms")
print("   📊 Figure 6: Detailed Correlation Matrix")
print("   📈 Figure 7: Strong Correlations Network")
print("\n🔍 Key Findings:")
print(f"   - Analyzed correlations between {len(correlation_df)} gene pairs")
print(f"   - Found {len(correlation_df[correlation_df['Abs_Correlation'] > 0.7])} strong correlations (|r| > 0.7)")
print(f"   - Top correlation: {correlation_df.iloc[0]['Gene1']}-{correlation_df.iloc[0]['Gene2']} (r = {correlation_df.iloc[0]['Correlation']:.3f})")
