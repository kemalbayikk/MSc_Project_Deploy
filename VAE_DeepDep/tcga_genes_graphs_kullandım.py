# import pandas as pd
# import matplotlib.pyplot as plt
# import os

# # Veriyi yükle
# data_path = 'Pytorch_codes/Variational_autoencoder/Models To Analyze/Split 2 VAE/tcga_predicted_data_best_model_vae_split_2.txt'  # Dosyanızın yolu
# data = pd.read_csv(data_path, sep='\t')

# cancer_type = "Lung_adenocarcinoma"
# # TCGA kodlarının hangi kanser türüne ait olduğunu gösteren dosyayı yükle
# cancer_type_path = f'Data/CancerTCGAMappings/{cancer_type}.txt'  # Breast invasive kanser tipi dosyası
# with open(cancer_type_path, 'r') as file:
#     breast_invasive_tcga_ids = file.read().splitlines()

# # Breast invasive kanser tipi ile ilgili kolonları filtrele
# filtered_data = data[['CRISPR_GENE'] + breast_invasive_tcga_ids]

# print(filtered_data)

# # Genlerin ortalama dağılımını hesapla
# filtered_data.set_index('CRISPR_GENE', inplace=True)
# mean_distribution = filtered_data.mean(axis=1)

# # En düşük 50 gene göre sırala
# lowest_50_genes = mean_distribution.nsmallest(20)

# # Grafiği oluştur ve kaydet
# output_dir = 'PytorchStaticSplits/DeepDepVAE/Analysis/TCGA Graphs'
# os.makedirs(output_dir, exist_ok=True)

# plt.figure(figsize=(20, 8))
# lowest_50_genes.plot(kind='bar')
# plt.title(f'The 20 Genes Lung Adenocarcinoma Tumors Are Most Dependent On', fontsize=28)
# plt.xlabel('Genes', fontsize=28)
# plt.ylabel('Mean Dependency Score', fontsize=28)
# plt.xticks(rotation=45, fontsize=18)
# plt.yticks(fontsize=14) 
# plt.tight_layout()
# plt.show()
# #plt.savefig(os.path.join(output_dir, f'lowest_50_genes_{cancer_type}_best.png'))
# plt.close()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Veriyi yükle
data_path = 'PytorchStaticSplits/DeepDepVAE/Analysis/tcga_predicted_data.txt'  # Dosyanızın yolu
data = pd.read_csv(data_path, sep='\t')

cancer_type = "Brain_Lower_Grade_Glioma"
# TCGA kodlarının hangi kanser türüne ait olduğunu gösteren dosyayı yükle
cancer_type_path = f'Data/CancerTCGAMappings/{cancer_type}.txt'  # Lung adenocarcinoma kanser tipi dosyası
with open(cancer_type_path, 'r') as file:
    lung_adenocarcinoma_tcga_ids = file.read().splitlines()

# Lung adenocarcinoma kanser tipi ile ilgili kolonları filtrele
filtered_data = data[['CRISPR_GENE'] + lung_adenocarcinoma_tcga_ids]

# Genlerin ortalama dağılımını hesapla ve en düşük 20 gene göre sırala
filtered_data.set_index('CRISPR_GENE', inplace=True)
mean_distribution = filtered_data.mean(axis=1)
lowest_20_genes = mean_distribution.nsmallest(20).index

# Sadece en düşük 20 geni içeren veriyi filtrele
lowest_20_genes_data = filtered_data.loc[lowest_20_genes]

# Veriyi long-formata dönüştür
lowest_20_genes_data = lowest_20_genes_data.reset_index().melt(id_vars='CRISPR_GENE', var_name='TCGA_ID', value_name='Dependency Score')

# Grafiği oluştur ve kaydet
output_dir = 'PytorchStaticSplits/DeepDepVAE/Analysis/TCGA Graphs'
os.makedirs(output_dir, exist_ok=True)

plt.figure(figsize=(20, 8))
sns.boxplot(x='CRISPR_GENE', y='Dependency Score', data=lowest_20_genes_data)
plt.title(f'The 20 Genes Brain Lower Grade Glioma Tumors Are Most Dependent On', fontsize=28)
plt.xlabel('Genes', fontsize=28)
plt.ylabel('Dependency Score', fontsize=28)
plt.xticks(rotation=45, fontsize=18)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()
#plt.savefig(os.path.join(output_dir, f'lowest_20_genes_boxplot_{cancer_type}_best.png'))
plt.close()

