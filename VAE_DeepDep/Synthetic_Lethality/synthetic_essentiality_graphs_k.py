# import matplotlib.pyplot as plt
# import pandas as pd

# gene = "BRCA1"

# # Load the Highest_SE_Scores_BRCA1.csv file
# se_scores_df = pd.read_csv(F'PytorchStaticSplits/SyntheticEssentiality/Highest_SE_Scores_{gene}_CCL.csv')

# # SE skorlarına göre sırala ve yalnızca en düşük 10 unique geni al
# df_unique = se_scores_df.drop_duplicates(subset='Gene')
# df_filtered = df_unique[df_unique['Gene'] != gene] 
# lowest_scores_df = df_filtered.sort_values(by='SE_Score').head(20)

# # Çizim
# plt.figure(figsize=(10, 6))
# plt.barh(f"{lowest_scores_df['Gene']}", lowest_scores_df['SE_Score'], color='skyblue')
# plt.xlabel('SE Score')
# plt.ylabel('Gene')
# plt.title('Top 20 Most Synthetic Essential Genes with BRCA1 Mutation (Lowest SE Scores)')
# plt.gca().invert_yaxis()  # En düşük skorlara göre yukarıdan aşağıya sıralama
# plt.show()

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

gene = "TP53"

# Load the Highest_SE_Scores_PTEN.csv file
se_scores_df = pd.read_csv(F'PytorchStaticSplits/SyntheticEssentiality/Highest_SE_Scores_{gene}_TCGA.csv')

# SE skorlarına göre sırala ve yalnızca en düşük 10 unique geni al
df_unique = se_scores_df.drop_duplicates(subset='Gene')
df_filtered = df_unique[df_unique['Gene'] != gene] 
lowest_scores_df = df_filtered.sort_values(by='SE_Score').head(10)

# En düşük skorlara sahip 10 geni al
lowest_genes = lowest_scores_df['Gene'].unique()

# Tüm skorlara göre genlerin plotu için filtreleme
plot_data = se_scores_df[se_scores_df['Gene'].isin(lowest_genes)]

gene_order = plot_data.groupby('Gene')['SE_Score'].mean().sort_values().index

# Çizim
plt.figure(figsize=(20, 8))
sns.boxplot(x='Gene', y='SE_Score', data=plot_data, order=gene_order, palette='tab10')
plt.ylabel('SL Score', fontsize=26)
plt.xlabel('Gene', fontsize=28)
plt.xticks(rotation=45, fontsize=24)
plt.yticks(fontsize=18) 
plt.title(f'Distribution of SL Scores for Top 10 Most Synthetic Lethal Genes with {gene} Mutation', fontsize=26)
#plt.gca().invert_xaxis()  # En düşük skorlara göre yukarıdan aşağıya sıralama
plt.tight_layout()
plt.show()
# # Her gen için SE skorlarının ortalamasını al, ancak kendisini hesaba katma
# avg_se_scores_df = df_filtered.groupby('Gene')['SE_Score'].mean().reset_index(name='Average_SE_Score')

# # En düşük 20 geni al
# lowest_scores_df = avg_se_scores_df.sort_values(by='Average_SE_Score').head(20)

# # Çizim
# plt.figure(figsize=(20, 8))
# sns.barplot(x='Average_SE_Score', y='Gene', data=lowest_scores_df, palette='tab10')
# plt.xlabel('Average SE Score')
# plt.ylabel('Gene')
# plt.title('Top 20 Most Synthetic Essential Genes with BRCA1 Mutation (Lowest Average SE Scores)')
# plt.gca().invert_yaxis()  # En düşük skorlara göre yukarıdan aşağıya sıralama
# plt.tight_layout()
#plt.show()