import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter

HERE = os.path.abspath(os.path.dirname(__file__))
DOWNLOAD_PATH = os.path.join(HERE, 'ArqFasta')
COM_CODONS = True


# ---------------------------------------------------------------------------- #
#                              LEITURA DOS DADOS                               #
# ---------------------------------------------------------------------------- #

def carregar_dados(download_path: str) -> list:
    identifier_pattern = re.compile(
        r"(?P<ID>[^\s]+)\s(?P<Host1>[A-Z])\/(?P<Host2>[^/]+(?:\s[^/]+)*)\/(?P<Location>[^\/]+)\/(?P<Data>[^/]+)\/(?P<Sequence_type>[^\s]+)\s"
    )

    data = []
    for file_name in os.listdir(download_path):
        if not file_name.endswith(".fa"):
            continue

        file_path = os.path.join(download_path, file_name)

        class_name = re.search(r"\((.*?)\)", file_name)
        class_name = class_name.group(1) if class_name else "Unknown"

        segment = re.search(r"_(.*?)_", file_name)
        segment = segment.group(1) if segment else "Unknown"

        hxnx = re.search(r"^(H\d+N\d+)", file_name)
        hxnx = hxnx.group(1) if hxnx else "Unknown"

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

    return data, identifier_pattern


# ---------------------------------------------------------------------------- #
#                            LIMPEZA E ESTRUTURAÇÃO                            #
# ---------------------------------------------------------------------------- #

def construir_dataframe(data: list, identifier_pattern) -> pd.DataFrame:
    data_cleaned = []
    for entry in data:
        match = identifier_pattern.match(entry[0])
        if match:
            data_cleaned.append([
                match.group("ID"),
                match.group("Host1"),
                match.group("Host2"),
                match.group("Location"),
                match.group("Data"),
                match.group("Sequence_type"),
                entry[1], entry[2], entry[3], entry[4]
            ])

    df = pd.DataFrame(data_cleaned, columns=[
        "ID", "Host1", "Host2", "Location", "Sequence_type", "Data",
        "Sequence", "Class", "Segment", "HxNx"
    ])

    df["Data"] = df["Data"].apply(lambda x: x[-4:] if len(x) > 4 else x)
    return df


def corrigir_colunas(df: pd.DataFrame) -> pd.DataFrame:
    df[["Data", "Sequence_type", "Location", "Host2"]] = df.apply(
        lambda row: pd.Series([
            row["Sequence_type"] if "/" in row["Data"] else row["Data"],
            row["Location"]      if "/" in row["Data"] else row["Sequence_type"],
            row["Host2"]         if "/" in row["Data"] else row["Location"],
            np.nan               if "/" in row["Data"] else row["Host2"]
        ]),
        axis=1
    )
    return df

def limpar_datas(df: pd.DataFrame) -> pd.DataFrame:
    df["Data"] = df["Data"].apply(lambda x: x[-4:] if len(x) > 4 else x)
    df["Data"] = df["Data"].str.replace(r"^([0-9]{2})$", r"19\1", regex=True)
    df = df[df["Data"].str.len() != 1]
    df = df[df["Data"].str.len() != 3]
    df["Data"] = df["Data"].apply(lambda x: x[1:] if len(x) == 5 else x)
    df["Data"] = df["Data"].apply(lambda x: x if x.isdigit() else None)
    df = df.dropna(subset=["Data"])
    return df

def filtrar_sequencias_validas(df: pd.DataFrame) -> pd.DataFrame:
    df["Not_AGCT"] = df["Sequence"].apply(
        lambda seq: np.nan if all(c in "AGCT" for c in seq)
        else sum(1 for c in seq if c not in "AGCT")
    )
    df = df[df["Not_AGCT"].isna()].copy()
    df["Nucleotides_count"] = df["Sequence"].apply(len)
    df["Codons_count"] = df["Sequence"].apply(lambda seq: len(seq) // 3)
    return df

# ---------------------------------------------------------------------------- #
#                             ANÁLISE DE CÓDONS                                #
# ---------------------------------------------------------------------------- #

def contar_codons(seq: str) -> dict:
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
    return dict(Counter(codons))

def gerar_tabela_codons(df: pd.DataFrame) -> pd.DataFrame:
    codon_dicts = df["Sequence"].apply(contar_codons)
    codon_df = pd.DataFrame.from_records(codon_dicts).fillna(0).astype(int)
    return codon_df.reset_index(drop=True)

def salvar_estatisticas_codons(codon_df: pd.DataFrame, segment: str, download_path: str):
    estatisticas = pd.DataFrame({
        "total": codon_df.sum(),
        "Media": codon_df.mean(),
        "DesvioPadrao": codon_df.std()
    })
    estatisticas.to_csv(os.path.join(download_path, f"estatisticas_codons_{segment}.csv"))

# ---------------------------------------------------------------------------- #
#                               VISUALIZAÇÕES                                  #
# ---------------------------------------------------------------------------- #

def plotar_distribuicao_por_data(df: pd.DataFrame):
    df_counts = df.groupby(["Data", "Segment"]).size().unstack(fill_value=0)
    df_counts.plot(kind="bar", stacked=False, color=["blue", "red"])
    plt.title("Distribuição de Sequências por Ano e Segmento")
    plt.xlabel("Ano")
    plt.ylabel("Contagem")
    plt.tight_layout()
    plt.show()

def plotar_boxplot_codons(codon_df: pd.DataFrame, segment: str, titulo: str):
    codon_df_sorted = codon_df[sorted(codon_df.columns)]

    mean_values = codon_df_sorted.mean()
    std_values  = codon_df_sorted.std()

    sns.boxplot(data=codon_df_sorted, color="lightgray", boxprops=dict(facecolor='gray'))

    for i, (mean, std) in enumerate(zip(mean_values, std_values)):
        plt.text(i, std,  '__', ha='center', va='center', color='blue', fontsize=10)
        plt.text(i, mean, '__', ha='center', va='center', color='red',  fontsize=10)

    media_patch   = mpatches.Patch(color='red',   label='Média')
    dp_patch      = mpatches.Patch(color='blue',  label='Desvio Padrão')
    mediana_patch = mpatches.Patch(color='black', label='Mediana')
    plt.legend(handles=[media_patch, dp_patch, mediana_patch], loc='upper right', bbox_to_anchor=(1, 1))

    plt.title(titulo)
    plt.xlabel("Códons")
    plt.ylabel("Contagem")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

def plotar_boxplot_simples(codon_df: pd.DataFrame, titulo: str):
    codon_df_sorted = codon_df[sorted(codon_df.columns)]
    sns.boxplot(data=codon_df_sorted, color="lightgray", boxprops=dict(facecolor='gray'))
    plt.title(titulo)
    plt.xlabel("Códons")
    plt.ylabel("Contagem")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

def plotar_comparacao_media(com_out: pd.DataFrame, sem_out: pd.DataFrame, segment: str):
    comparacao = pd.DataFrame({
        "sem_out": sem_out.mean(),
        "com_out": com_out.mean()
    })
    comparacao.plot(kind="bar", color=["gray", "black"], alpha=0.8)
    plt.xlabel("Códons")
    plt.ylabel("Contagem Média")
    plt.title(f"Comparação de Média de Códons com e sem Outliers ({segment})")
    plt.xticks(rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.show()

# ---------------------------------------------------------------------------- #
#                           DETECÇÃO DE OUTLIERS                               #
# ---------------------------------------------------------------------------- #

def detectar_outliers(df: pd.DataFrame, codon_df: pd.DataFrame) -> tuple:
    Q1  = codon_df.quantile(0.25)
    Q3  = codon_df.quantile(0.75)
    IQR = Q3 - Q1
    outlier_mask = (codon_df < (Q1 - 1.5 * IQR)) | (codon_df > (Q3 + 1.5 * IQR))

    outlier_df = codon_df[outlier_mask]
    tabela_raw = pd.concat([df["ID"].reset_index(drop=True), outlier_df.reset_index(drop=True)], axis=1)

    # Formata: "valor(codon), valor(codon), ..."
    def formatar_linha(row):
        return [f"{val}({col})" for col, val in row.items() if pd.notna(val)]

    linhas = tabela_raw.apply(formatar_linha, axis=1)
    Tabela_De_Outliers = pd.DataFrame({"Valores": [", ".join(l) for l in linhas]})
    Tabela_De_Outliers = Tabela_De_Outliers[Tabela_De_Outliers["Valores"].str.count(r"\(") >= 2]

    df_out_temp = pd.concat([df["ID"].reset_index(drop=True), outlier_df.reset_index(drop=True)], axis=1)
    df_out_temp["tem_out"] = (~df_out_temp.drop(columns="ID").isna().all(axis=1)).astype(int)

    df_com_out = df_out_temp[df_out_temp["tem_out"] == 1].drop(columns="tem_out")
    df_sem_out = df_out_temp[df_out_temp["tem_out"] == 0].drop(columns="tem_out")

    return Tabela_De_Outliers, df_com_out, df_sem_out, outlier_mask

def codons_por_grupo(df: pd.DataFrame, codon_df: pd.DataFrame, ids: pd.Series, segment: str) -> pd.DataFrame:
    df_filtrado = df[df["ID"].isin(ids)].reset_index(drop=True)
    resultado = pd.concat([df_filtrado[["ID", "Segment"]], codon_df.reset_index(drop=True)], axis=1)
    resultado = resultado.dropna(subset=["ID"])
    return resultado[resultado["Segment"] == segment]

# ---------------------------------------------------------------------------- #
#                           PRINTS E RELATÓRIOS                                #
# ---------------------------------------------------------------------------- #

def printar_secao(titulo: str):
    borda = "-" * 78
    print(f"\n{borda}")
    print(f"{titulo:^78}")
    print(borda)

def relatorio_inicios_fins(df: pd.DataFrame):
    printar_secao("Inícios e fins específicos")
    print(f"Sequências que começam com 'ATG': {df['Sequence'].str.startswith('ATG').sum()}")
    print(f"Sequências que terminam com 'TAA': {df['Sequence'].str.endswith('TAA').sum()}")
    print(f"Sequências que terminam com 'TAG': {df['Sequence'].str.endswith('TAG').sum()}")
    print(f"Sequências que terminam com 'TGA': {df['Sequence'].str.endswith('TGA').sum()}")
    total_stop = (df["Sequence"].str.endswith("TAA") |
                  df["Sequence"].str.endswith("TAG") |
                  df["Sequence"].str.endswith("TGA")).sum()
    print(f"Total terminadas com TAG, TGA ou TAA: {total_stop}")

def main():
    data, identifier_pattern = carregar_dados(DOWNLOAD_PATH)

    df = construir_dataframe(data, identifier_pattern)
    df.to_csv(os.path.join(DOWNLOAD_PATH, "csvcerto2.csv"), index=False)

    df = corrigir_colunas(df)
    df = limpar_datas(df)

    plotar_distribuicao_por_data(df)

    df = filtrar_sequencias_validas(df)

    printar_secao("DataFrame básico")
    print(df)

    if not COM_CODONS:
        return

    codon_df = gerar_tabela_codons(df)
    df_codons = pd.concat([df.reset_index(drop=True), codon_df], axis=1)

    codon_usage_table = pd.concat([df["ID"].reset_index(drop=True), codon_df], axis=1)
    codon_usage_table.to_csv(os.path.join(DOWNLOAD_PATH, "codon_usage_table.csv"), index=False)

    for seg in ["NA", "HA"]:
        subset = df_codons[df_codons["Segment"] == seg].drop(columns=df.columns, errors='ignore')
        salvar_estatisticas_codons(subset, seg, DOWNLOAD_PATH)
        plotar_boxplot_codons(
            codon_df=subset,
            segment=seg,
            titulo=f"Box Plot com Média e Desvio Padrão dos Códons — Influenza A ({seg})"
        )

  
    relatorio_inicios_fins(df)

    Tabela_De_Outliers, df_com_out, df_sem_out, outlier_mask = detectar_outliers(df, codon_df)

    Tabela_De_Outliers.to_csv(os.path.join(DOWNLOAD_PATH, "Tabela_De_Outliers.csv"), index=False)
    df_com_out.to_csv(os.path.join(DOWNLOAD_PATH, "df_com_out.csv"), index=False)
    df_sem_out.to_csv(os.path.join(DOWNLOAD_PATH, "df_sem_out.csv"), index=False)

    printar_secao("Tabela de Outliers")
    print(Tabela_De_Outliers)

    printar_secao("Sequências com outliers")
    print(df_com_out)
    df_com_out.info()

    printar_secao("Sequências sem outliers")
    print(df_sem_out)
    df_sem_out.info()

    for seg in ["HA", "NA"]:
        com_out = codons_por_grupo(df, codon_df, df_com_out["ID"], seg)
        sem_out = codons_por_grupo(df, codon_df, df_sem_out["ID"], seg)

        printar_secao(f"Boxplot — sequências COM outliers ({seg})")
        plotar_boxplot_simples(
            com_out.drop(columns=["ID", "Segment"], errors='ignore'),
            titulo=f"Box Plot dos Códons — Influenza A ({seg}) com Outliers"
        )
        print(com_out)
        com_out.info()

        printar_secao(f"Boxplot — sequências SEM outliers ({seg})")
        plotar_boxplot_simples(
            sem_out.drop(columns=["ID", "Segment"], errors='ignore'),
            titulo=f"Box Plot dos Códons — Influenza A ({seg}) sem Outliers"
        )
        print(sem_out)

        printar_secao(f"Comparação com e sem outliers ({seg})")
        plotar_comparacao_media(
            com_out=com_out.drop(columns=["ID", "Segment"], errors='ignore'),
            sem_out=sem_out.drop(columns=["ID", "Segment"], errors='ignore'),
            segment=seg
        )

    printar_secao("Gerando arquivo FASTA com sequências com outliers")
    sequencias_com_out = df[df["ID"].isin(df_com_out["ID"])]
    output_path = os.path.join(DOWNLOAD_PATH, "filtrado.fa")
    with open(output_path, "w") as f:
        for _, row in sequencias_com_out.iterrows():
            f.write(f">{row['ID']}\n{row['Sequence']}\n")
    print(f"Arquivo FASTA salvo em: {output_path}")

if __name__ == "__main__":
    main()