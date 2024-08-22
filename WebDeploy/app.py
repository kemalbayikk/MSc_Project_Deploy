from flask import Flask, render_template, request, jsonify
import matplotlib
matplotlib.use('Agg')  # GUI olmadan matplotlib kullanımı için
import io
import base64
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import numpy as np
import pandas as pd
from V_DeepDep import DeepDEP
import threading
import time

app = Flask(__name__)

# Global değişken ilerleme durumunu ve sonucu tutar
progress_data = {
    "progress": 0,
    "table_html_lowest": None,
    "table_html": None,
    "error_message": None  # Hata mesajını tutmak için
}

def long_running_task(selected_dataset, gene):
    global progress_data
    print(f"Thread started with gene: {gene}")
    progress_data['progress'] = 0  # İlerlemeyi sıfırla
    progress_data['error_message'] = None  # Hata mesajını sıfırla

    fprint_dim = 3115
    dense_layer_dim = 100

    # Dataset seçimine göre veri yükleme
    if selected_dataset == "CCL":
        data_mut, _, sample_names_mut, gene_names_mut = load_data("Data/ccl_mut_data_paired_with_tcga.txt")
    else:
        data_mut, _, sample_names_mut, gene_names_mut = load_data("Data/tcga_mut_data_paired_with_ccl.txt")

    gene = gene.upper()
    
    # Gen listesinde aranan genin olup olmadığını kontrol et
    if gene not in gene_names_mut:
        progress_data['error_message'] = f"Error: The gene '{gene}' was not found in the dataset."
        print(progress_data['error_message'])
        return

    data_fprint_1298DepOIs, _, gene_names_fprint, _ = load_data("Data/crispr_gene_fingerprint_cgp.txt")
    dims_mut = 4539
    model = DeepDEP(dims_mut, fprint_dim, dense_layer_dim)
    model.load_state_dict(torch.load(f"VAE_DeepDep/DeepDepMutationOnly/Results/Split2/PredictionNetworkModels/VAE_Prediction_Network_Split_2_Only_Mutation.pth"))
    model.eval()

    gene_index = gene_names_mut.index(gene)
    data_pred = np.zeros((len(sample_names_mut), data_fprint_1298DepOIs.shape[0]))
    max_sl_scores = []

    for z in range(len(sample_names_mut)):
        progress_data['progress'] = int((z + 1) / len(sample_names_mut) * 99)  # İlerleme yüzdesini güncelle
        data_mut_batch = np.zeros((data_fprint_1298DepOIs.shape[0], dims_mut), dtype='float32')
        data_mut_batch[:, :] = data_mut[z, :]
        data_mut_batch = torch.tensor(data_mut_batch, dtype=torch.float32)
        data_fprint_batch = torch.tensor(data_fprint_1298DepOIs, dtype=torch.float32)

        with torch.no_grad():
            data_mut_batch[:, gene_index] = 1.0
            output_mut = model(data_mut_batch, data_fprint_batch).cpu().numpy()
            data_mut_batch[:, gene_index] = 0.0
            output_wt = model(data_mut_batch, data_fprint_batch).cpu().numpy()

        sl_scores = output_mut - output_wt
        data_pred[z] = np.transpose(sl_scores)
        min_se_index = np.argmin(sl_scores)
        max_sl_scores.append((sample_names_mut[z], gene_names_fprint[min_se_index], sl_scores[min_se_index][0]))

    print(f"Task completed with gene: {gene}")

    # Grafik oluşturma
    max_sl_scores_sorted = sorted(max_sl_scores, key=lambda x: x[2])
    max_sl_scores_df = pd.DataFrame(max_sl_scores_sorted, columns=['Tumor', 'Gene', 'Synthetic Lethality Score'])
    df_unique = max_sl_scores_df.drop_duplicates(subset='Gene')
    df_filtered = df_unique[df_unique['Gene'] != gene]
    lowest_scores_df = df_filtered.sort_values(by='Synthetic Lethality Score').head(10)

    # lowest_genes = lowest_scores_df['Gene'].unique()
    # plot_data = max_sl_scores_df[max_sl_scores_df['Gene'].isin(lowest_genes)]
    # gene_order = plot_data.groupby('Gene')['SL_Score'].mean().sort_values().index

    # plt.figure(figsize=(16, 10))
    # sns.boxplot(x='Gene', y='SL_Score', data=plot_data, order=gene_order, palette='tab10')
    # plt.ylabel('SL Score', fontsize=18)
    # plt.xlabel('Gene', fontsize=18)
    # plt.xticks(rotation=45, fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.title(f'Distribution of SL Scores for Top 10 Most Synthetic Lethal Genes with {gene} Mutation', fontsize=18)

    # img = io.BytesIO()
    # plt.savefig(img, format='png')
    # img.seek(0)
    # progress_data['plot_url'] = base64.b64encode(img.getvalue()).decode()  # Grafiği global değişkene kaydet

    # DataFrame oluşturma (data_pred verilerini tabloya aktarma)
    data_pred_df = pd.DataFrame(data_pred, index=sample_names_mut, columns=gene_names_fprint)

    data_pred_df = data_pred_df.drop(gene, axis=1, errors='ignore')  # Sütun olarak geni çıkar
    data_pred_df = data_pred_df.drop(gene, axis=0, errors='ignore')  # Satır olarak geni çıkar

    data_pred_df = data_pred_df.transpose()

    # DataFrame'i HTML formatına çevirme
    progress_data['table_html'] = data_pred_df.to_html(classes='table table-striped', table_id='synthetic-table')
    # Tabloyu HTML formatına çevir ve global değişkene kaydet
    
    max_sl_scores_df = max_sl_scores_df[max_sl_scores_df['Gene'] != gene]
    progress_data['table_html_lowest'] = max_sl_scores_df.sort_values(by='Synthetic Lethality Score').to_html(classes='table table-striped', table_id='synthetic-table')
    progress_data['progress'] = 100  # İşlem tamamlandı

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/progress', methods=['GET'])
def progress_status():
    global progress_data
    return jsonify({
        'progress': progress_data['progress'],
        'error_message': progress_data['error_message']  # Hata mesajını da döndürüyoruz
    })

@app.route('/result', methods=['GET'])
def result():
    global progress_data
    table_html_lowest = progress_data['table_html_lowest']
    table_html = progress_data['table_html']
    error_message = progress_data['error_message']

    if error_message:
        return f'<div class="alert alert-danger" role="alert">{error_message}</div>'
    elif table_html:
        return f'''
            <div class="mt-5">
                <h3>All Synthetic Lethality Scores</h3>
                {table_html}
            </div>
            <div class="mt-5">
                <h3>Lowest Synthetic Lethality Scores</h3>
                {table_html_lowest}
            </div>
        '''
    else:
        return "Prediction still in progress, please check again later."

@app.route('/predict', methods=['POST'])
def predict():
    selected_dataset = request.form['dataset']
    selected_gene = request.form['gene']
    thread = threading.Thread(target=long_running_task, args=(selected_dataset, selected_gene))
    thread.start()  # İşlemi arka planda başlat
    print(f"Prediction thread started for gene: {selected_gene}")

    # İşlem tamamlanana kadar bekle
    # thread.join()

    return render_template('index.html')
    
# İkinci tab için dosya okuma işlemi (Dependency Predictions)
@app.route('/dependency_predictions', methods=['POST'])
def dependency_predictions():
    selected_dataset = request.form['dataset']
    gene = request.form['gene'].upper()
    # CCL dataset seçildiyse ccl_predictions.txt dosyasını göster
    if selected_dataset == "CCL":
        file_path = "VAE_DeepDep/Predictions/ccl_predicted_data.txt"  # Dosya yolunu düzenleyin
    elif selected_dataset == "TCGA":
        file_path = "VAE_DeepDep/Predictions/tcga_predicted_data.txt"  # TCGA dosyasının yolu
    else:
        return "Dataset not supported yet."
    
    try:
        # Dosyayı yükle
        df = pd.read_csv(file_path, delimiter='\t')
        
        # Girilen gene göre filtrele
        filtered_df = df[df['CRISPR_GENE'].str.upper() == gene]
        
        # Eğer sonuç yoksa kullanıcıya bilgilendirme mesajı döndür
        if filtered_df.empty:
            return f"No results found for gene: {gene}"
        
        # Sütun başlıklarını kaydedelim
        column_headers = list(filtered_df.columns)

        # Tabloyu transpoze et (satır ve sütunları yer değiştir)
        transposed_df = filtered_df.transpose()

        # Transpoze edilen tabloya başlık ekleyelim
        transposed_df.columns = [gene]

        # Başlıkları da manuel olarak ilk sütun olarak ekleyelim
        transposed_df.insert(0, 'Tumors', column_headers)

        transposed_df = transposed_df.drop(transposed_df.index[0])
        
        # Tabloyu HTML formatına çevir
        table_html = transposed_df.to_html(classes='table table-striped', index=False, table_id='dependency-table')
        return table_html  # HTML olarak tabloyu döndür
    except Exception as e:
        return f"Error loading file: {str(e)}"


def load_data(filename):
    data = []
    gene_names = []
    data_labels = []
    with open(filename) as file:
        lines = file.readlines()
        sample_names = lines[0].replace('\n', '').split('\t')[1:]
        for line in lines[1:]:
            values = line.replace('\n', '').split('\t')
            gene = values[0].upper()
            gene_names.append(gene)
            data.append(values[1:])
    data = np.array(data, dtype='float32')
    data = np.transpose(data)

    return data, data_labels, sample_names, gene_names

if __name__ == '__main__':
    app.run(port=3030, debug=True)
