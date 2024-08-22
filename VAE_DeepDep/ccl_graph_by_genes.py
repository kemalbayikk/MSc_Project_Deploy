import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Verileri yükleyelim
csv_data_path = 'Data/data_dep_updated_2.csv'  # CSV dosyası için düzgün yolu girin
txt_data_path = 'Pytorch_codes/Variational_autoencoder/Models To Analyze/Split 2 VAE/ccl_predicted_data_best_model_vae_split_2.txt'  # TXT dosyası için düzgün yolu girin

csv_data = pd.read_csv(csv_data_path)
txt_data = pd.read_csv(txt_data_path, delimiter='\t')  # Eğer dosya tab ile ayrılmışsa delimiter='\t' kullanın

# CRISPR_GENE sütununu koruyarak verileri geniş formatından uzun formata dönüştür
csv_melted = csv_data.melt(id_vars=['CRISPR_GENE'], var_name='CCL', value_name='Dependency Score').assign(Source='Original Data')
txt_melted = txt_data.melt(id_vars=['CRISPR_GENE'], var_name='CCL', value_name='Dependency Score').assign(Source='Predicted Data')

# Her iki veri seti için genlerin ortalama etkisini hesapla
csv_average_impact = csv_melted.groupby('CRISPR_GENE')['Dependency Score'].mean().reset_index()

# En etkili 30 geni seç
top_30_genes_csv = csv_average_impact.nsmallest(10, 'Dependency Score')['CRISPR_GENE']

# Veri setlerini birleştir
filtered_csv_data = csv_melted[csv_melted['CRISPR_GENE'].isin(top_30_genes_csv)]
filtered_txt_data = txt_melted[txt_melted['CRISPR_GENE'].isin(top_30_genes_csv)]

combined_data = pd.concat([filtered_csv_data, filtered_txt_data])

gene_order = combined_data.groupby('CRISPR_GENE')['Dependency Score'].mean().sort_values().index

# # Violin plot çiz
plt.figure(figsize=(20, 10))
sns.violinplot(x='CRISPR_GENE', y='Dependency Score', hue='Source', data=combined_data, split=True, palette=['blue', 'orange'], order=gene_order)
plt.title('Comparison of Dependency Scores for Top 10 Most Effective Genes Across All CCLs', fontsize=28)
plt.xticks(rotation=45, fontsize=28)  # Daha fazla gen ismini görmek için etiketleri döndür
plt.yticks(fontsize=14) 
plt.xlabel('Gene', fontsize=28)
plt.ylabel('Dependency Score', fontsize=28)
plt.legend(title='Data Source', loc='lower right', fontsize=18, title_fontsize=18)

# Violin plot çiz
# plt.figure(figsize=(20, 8))
# sns.boxplot(x='CRISPR_GENE', y='Dependency Score', hue='Source', data=combined_data, palette=['blue', 'orange'], order=gene_order)
# plt.title('Comparison of Dependency Scores for Top 10 Most Effective Genes Across All CCLs', fontsize=28)
# plt.xticks(rotation=45, fontsize=28)  # Daha fazla gen ismini görmek için etiketleri döndür
# plt.yticks(fontsize=14) 
# plt.xlabel('Gene', fontsize=28)
# plt.ylabel('Dependency Score', fontsize=28)
# plt.legend(title='Data Source', loc='lower right', fontsize=18, title_fontsize=18)

# Grafiği göster
plt.tight_layout()  # Layout'u düzenle
plt.show()