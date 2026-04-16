import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

df_analise_N = pd.read_pickle("df_analise_N1.pkl") #pkl passou as colunas com uma consistencia maior, sem erros



#----------------------------------------------------------------------#
#                  IDENTIFICAR CARACTERES NÃO ACGT                     #
#----------------------------------------------------------------------#
caracteres_usados = set("".join(df_analise_N["Sequence"]))
caracteres_diferentes = sorted(caracteres_usados - set("ACGT"))

print("Caracteres diferentes de ACGT encontrados:", caracteres_diferentes)

#----------------------------------------------------------------------#
#                   GERAÇÃO DE GRÁFICOS POR CARACTERE                  #
#----------------------------------------------------------------------#
for char in caracteres_diferentes:
    coluna_pct = f"Perc_{char}" 
    df_analise_N[coluna_pct] = df_analise_N["Sequence"].apply(lambda seq: seq.count(char) / len(seq))
    
    # Limites de 50% a 1%(valor maximo[50], valor minimo[0.9], passo[-1])
    limites = np.arange(50, 0.9, -1)
    resultados = []

    for limite in limites:
        limite_decimal = limite / 100
        num_removidas = (df_analise_N[coluna_pct] > limite_decimal).sum()
        resultados.append({"Limite_%": limite, "Removidas": num_removidas})

    df_resultados = pd.DataFrame(resultados)
    removidas_com_1 = (df_analise_N["Sequence"].str.contains(char)).sum()

   

   
plt.figure(figsize=(10, 5))
plt.plot(df_resultados["Limite_%"], df_resultados["Removidas"], marker='o', color='blue', label=f"Remoção por limite % de '{char}'")
plt.scatter(0, removidas_com_1, color='red', s=100, zorder=5, label=f"Remoção se houver ≥1 '{char}'")
plt.text(1, removidas_com_1, f"{removidas_com_1} seqs", color='red', va='bottom')
plt.title(f"Remoção de Sequências com '{char}'")
plt.xlabel("Limite máximo (%) do caractere na sequência")
plt.ylabel("Número de sequências removidas")
#plt.gca().invert_xaxis()
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

df_analise_N["Quantidade_caracteres"] = df_analise_N["Sequence"].str.len()
df_analise_N["Quantidade_caracteres_estranhos"] = df_analise_N["Sequence"].apply(lambda x: sum(1 for c in x.upper() if c not in "ACGT"))

resultados_gerais = []

for limite in limites: 
        limite_real = limite / 100
        num_removidas = (
            df_analise_N["Quantidade_caracteres_estranhos"] > 
            (df_analise_N["Quantidade_caracteres"] * limite_real)).sum()
        resultados_gerais.append({"Limite_%": limite, "Removidas": num_removidas})

df_resultados_gerais = pd.DataFrame(resultados_gerais)
removidas_com_1 = len(df_analise_N)

plt.figure(figsize=(10, 5))
plt.plot(df_resultados_gerais["Limite_%"], df_resultados_gerais["Removidas"],marker='o', color='purple', label="Remoção por % total de caracteres estranhos")
plt.scatter(0, removidas_com_1, color='red', s=100, zorder=5, label="Remoção se houver ≥1 caractere estranho")
plt.text(1, removidas_com_1, f"{removidas_com_1} seqs", color='red', va='bottom')

plt.title("Remoção por Percentual Total de Caracteres Estranhos")
plt.xlabel("Limite máximo (%) de caracteres não ACGT")
plt.ylabel("Número de sequências removidas")
#plt.gca().invert_xaxis()
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()      
print(df_analise_N[["Sequence", "Quantidade_caracteres", "Quantidade_caracteres_estranhos"]].head(10))

plt.figure(figsize=(10, 5))
plt.plot(df_resultados_gerais["Limite_%"], df_resultados_gerais["Removidas"],marker='o', color='purple', label="Remoção por % total de caracteres estranhos")
plt.title("Remoção por Percentual Total de Caracteres Estranhos")
plt.xlabel("Limite máximo (%) de caracteres não ACGT")
plt.ylabel("Número de sequências removidas")
    #plt.gca().invert_xaxis()
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()      
print(df_analise_N[["Sequence", "Quantidade_caracteres", "Quantidade_caracteres_estranhos"]].head(10))

#----------------------------------------------------------------------#
#                        FILTRANDO SEQUÊNCIAS                          #
#----------------------------------------------------------------------#
df_analise_N = df_analise_N[df_analise_N["Sequence"].str.contains("[^ACGT]")].copy()


def caracteres_nao_acgt(seq):
    return ''.join(sorted(set(seq) - set("ACGT")))
def quantidade_caracteres(seq):
    return len(seq)
def quantidade_caracteres_estranhos(seq):
    return sum(1 for c in seq if c not in "ACGT")
df_analise_N["Caracteres_estranhos"] = df_analise_N["Sequence"].apply(caracteres_nao_acgt)
df_analise_N["Quantidade_caracteres"] = df_analise_N["Sequence"].apply(quantidade_caracteres)
df_analise_N["Quantidade_caracteres_estranhos"] = df_analise_N["Sequence"].apply(quantidade_caracteres_estranhos)

print("filtrar sequências por quantidade percentual de caracteres estranhos?")
resposta_seq = input().strip().lower()
while (resposta_seq == 's'):
    print("qual o limite máximo percentual de caracteres estranhos?")
    limite_percentual = float(input().strip())

    filtro = df_analise_N["Quantidade_caracteres_estranhos"] < (df_analise_N["Quantidade_caracteres"] * limite_percentual / 100)
    print(df_analise_N.loc[filtro, ["Sequence", "Caracteres_estranhos", "Quantidade_caracteres", "Quantidade_caracteres_estranhos"]])
    

    
    df_ha_na = df_analise_N.loc[filtro & df_analise_N["Segment"].isin(["HA", "NA"])]

    
    media_tamanhos = df_ha_na.groupby("Segment")["Quantidade_caracteres"].mean()
    print("Média de tamanhos por segmento (HA vs NA):")
    print(media_tamanhos)
    #----------------------------------------------------------------------#
    #             MEDIA DE TAMANHO DAS SEQUÊNCIAS POR SEGMENTO             # 
    #----------------------------------------------------------------------#
    plt.figure(figsize=(8, 5))
    media_tamanhos.plot(kind="bar", color=["skyblue", "lightcoral"])
    plt.title(f"Tamanho médio por segmento (HA vs NA) — Limite: {limite_percentual}%")
    plt.ylabel("Tamanho médio")
    plt.xlabel("Segmento")
    plt.xticks(rotation=0)
    plt.grid(axis='y')
    plt.tight_layout()
    plt.show()
    #----------------------------------------------------------------------#
    #                  BOX PLOT DE TAMANHO DAS SEQUÊNCIAS                  # 
    #----------------------------------------------------------------------#
    df_filtrado = df_analise_N.loc[filtro]  
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_filtrado, x="Segment", y="Quantidade_caracteres", palette="Set2")
    plt.title(f"Boxplot: Tamanho das Sequências por Segmento — Limite: {limite_percentual}%")
    plt.xlabel("Segmento")
    plt.ylabel("Tamanho da Sequência")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.show()
    print("continuar filtrando sequências por quantidade percentual de caracteres estranhos?")
    resposta_seq = input().strip().lower()
print(df_analise_N[["Sequence", "Caracteres_estranhos", "Quantidade_caracteres", "Quantidade_caracteres_estranhos"]])






