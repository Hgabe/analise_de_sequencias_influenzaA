import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import matplotlib.patches as mpatches

download_path = r"D:\codigopython\ArqFasta"
comcodons = True
data = []

# Regex ajustado para extrair corretamente Host1 e Host2
identifier_pattern = re.compile(
    r"(?P<ID>[^\s]+)\s(?P<Host1>[A-Z])\/(?P<Host2>[^/]+(?:\s[^/]+)*)\/(?P<Location>[^\/]+)\/(?P<Data>[^/]+)\/(?P<Sequence_type>[^\s]+)\s"
)

for file_name in os.listdir(download_path):
    if file_name.endswith(".fa"):
        file_path = os.path.join(download_path, file_name)

        # Captura do nome da classe assim: (x)
        class_match = re.search(r"\((.*?)\)", file_name)
        class_name = class_match.group(1) if class_match else "Unknown"

        # Captura do nome do segmento assim: _x_
        segment_match = re.search(r"_(.*?)_", file_name)
        segment = segment_match.group(1) if segment_match else "Unknown"

        # Captura o HxNx no nome
        hxnx_match = re.search(r"^(H\d+N\d+)", file_name)
        hxnx = hxnx_match.group(1) if hxnx_match else "Unknown"

        # Abre o arquivo fasta
        with open(file_path, "r") as fasta_file:
            identifier = ""
            sequence = ""

            for line in fasta_file:
                line = line.strip()
                if line.startswith(">"):
                    if identifier and sequence:
                        data.append([identifier, sequence, class_name, segment, hxnx])
                    identifier = line[1:]
                    sequence = ""
                else:
                    sequence += line

            if identifier and sequence:
                data.append([identifier, sequence, class_name, segment, hxnx])

# Aplica regex para extrair colunas
data_cleaned = []
for entry in data:
    match = identifier_pattern.match(entry[0])
    if match:
        ID = match.group("ID")
        Host1 = match.group("Host1")
        Host2 = match.group("Host2")
        Location = match.group("Location")
        Data = match.group("Data")
        Sequence_type = match.group("Sequence_type")

        data_cleaned.append([
            ID, Host1, Host2, Location, Data, Sequence_type,
            entry[1], entry[2], entry[3], entry[4]
        ])

# Cria DataFrame
df = pd.DataFrame(data_cleaned, columns=[
    "ID", "Host1", "Host2", "Location", "Sequence_type", "Data",
    "Sequence", "Class", "Segment", "HxNx"
])

df["Data"] = df["Data"].apply(lambda x: x[-4:] if len(x) > 4 else x)


# Salva o CSV
csv_path = os.path.join(download_path, "csvcerto2.csv")
df.to_csv(csv_path, index=False)
print(f"Arquivo CSV salvo em: {csv_path}")
# Corrige colunas embaralhadas por erros de estrutura no identificador
df[["Data", "Sequence_type", "Location", "Host2"]] = df.apply(
    lambda row: pd.Series([
        row["Sequence_type"] if "/" in row["Data"] else row["Data"],
        row["Location"] if "/" in row["Data"] else row["Sequence_type"],
        row["Host2"] if "/" in row["Data"] else row["Location"],
        np.nan if "/" in row["Data"] else row["Host2"]
    ]),
    axis=1
)


# Correções de datas e remoção de registros inválidos
df["Data"] = df["Data"].apply(lambda x: x[-4:] if len(x) > 4 else x)
df["Data"] = df["Data"].str.replace(r"^([0-9]{2})$", r"19\1", regex=True)
df = df[df["Data"].str.len() != 1]
df = df[df["Data"].str.len() != 3]
df["Data"] = df["Data"].apply(lambda x: x[1:] if len(x) == 5 else x)
df["Data"] = df["Data"].apply(lambda x: x if x.isdigit() else None)
df = df.dropna(subset=["Data"])


# Geração do gráfico
df_counts = df.groupby(["Data", "Segment"]).size().unstack(fill_value=0)
df_counts.plot(kind="bar", stacked=False, color=["blue", "red"])
plt.show()

# Se a sequencia so tiver caracteres validos apresenta NaN 
df["Not_AGCT"] = df["Sequence"].apply(lambda seq: np.nan if all(c in "AGCT" for c in seq) else sum(1 for c in seq if c not in "AGCT"))
df = df[df["Not_AGCT"].isna()]
df["Nucleotides_count"] = df["Sequence"].apply(lambda seq: len (seq))
df["Codons_count"] = df["Sequence"].apply(lambda seq: len (seq)//3)
if comcodons:
    def contar_codons(seq):
        codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
        return dict(Counter(codons))

    # Aplica a função para cada linha da coluna "Sequence"
    codon_dicts = df["Sequence"].apply(contar_codons)

    # Transforma a lista de dicionários em colunas
    codon_df = pd.DataFrame.from_records(codon_dicts).fillna(0).astype(int)
    df_codons = pd.concat([df.reset_index(drop=True), codon_df.reset_index(drop=True)], axis=1)
    codon_dfs = pd.concat([df["Segment"].reset_index(drop=True), codon_df.reset_index(drop=True)], axis=1)
    
    def geraboxplot():
        estatisticas_codons = pd.DataFrame({
            "total": codon_df.sum(),
            "Media": codon_df.mean(),
            "DesvioPadrao": codon_dfs.std()
        })
        estatisticas_codons.to_csv(os.path.join(download_path, f"estatisticas_codons_{segment}.csv"))

        codon_df_sorted = codon_dfs[sorted(codon_df.columns)]
        sns.boxplot(data=codon_df_sorted, color="lightgray", boxprops=dict(facecolor='gray'))
        
        mean_values = codon_df_sorted.mean()
        std_values = codon_df_sorted.std()

        for i, (mean, std) in enumerate(zip(mean_values, std_values)):
            plt.text(i, std, '__', ha='center', va='center', color='blue', fontsize=10)
            plt.text(i, mean, '__', ha='center', va='center', color='red', fontsize=10)

        plt.title("Box Plot com Média e Desvio Padrão dos Códons de Sequências do Vírus Influenza A Para o Gene " + segment)
        plt.xlabel("Códons")
        plt.ylabel("Contagem")
        plt.xticks(rotation=90)

        media_patch = mpatches.Patch(color='red', label='Média')
        dp_patch = mpatches.Patch(color='blue', label='Desvio Padrão')
        mediana_patch = mpatches.Patch(color='black', label='mediana')
        plt.legend(handles=[media_patch, dp_patch, mediana_patch], loc='upper right', bbox_to_anchor=(1, 1))

        plt.tight_layout()
        plt.show()

    # Filtra e plota apenas para segmentos específicos
    segment = "NA"
    codon_dfs = df_codons[df_codons["Segment"] == segment]
    geraboxplot()

    segment = "HA"
    codon_dfs = df_codons[df_codons["Segment"] == segment]
    geraboxplot()

    df.info()
    
    
    # Junta apenas a coluna "ID" com as contagens dos códons
    codon_usage_table = pd.concat([df["ID"].reset_index(drop=True), codon_df.reset_index(drop=True)], axis=1)

    # Define o caminho para salvar o CSV
    codon_csv_path = os.path.join(download_path, "codon_usage_table.csv")

    # Salva em CSV
    codon_usage_table.to_csv(codon_csv_path, index=False)
    std_series = df_codons.drop(columns=df.columns).std()

    
    print(f"Arquivo com ID e contagens de códons salvo em: {codon_csv_path}")

print(df)
df.info()

count_atg = df["Sequence"].str.startswith("ATG").sum()
print(f"Número de sequências que começam com 'ATG': {count_atg}")
count_taa = df["Sequence"].str.endswith("TAA").sum()
print(f"Número de sequências que terminam com 'TAA': {count_taa}")
count_tag = df["Sequence"].str.endswith("TAG").sum()
print(f"Número de sequências que terminam com 'TAG': {count_tag}")
count_tga = df["Sequence"].str.endswith("TGA").sum()
print(f"Número de sequências que terminam com 'TGA': {count_tga}")
print(f"Número total de sequêncais terminadas com TAG,TGA,TAA: {count_taa + count_tag + count_tga}")

