## Descrição

Este projeto realiza uma análise bioinformática completa das sequências dos genes HA (Hemaglutinina) e NA (Neuraminidase) do vírus Influenza A, abrangendo hospedeiros aviário, humano e suíno.

O pipeline é dividido em quatro etapas principais:

1. **Download automatizado** das sequências do NCBI Influenza Virus Database
2. **Pré-processamento** e filtragem das sequências
3. **Análise exploratória** de caracteres inválidos
4. **Análise de uso de códons** e detecção de outliers

---

## Estrutura do Projeto

```
influenza-codon-analysis/
│
├── download.py              # Download automatizado via Selenium (NCBI Influenza DB)
├── preprocess.py            # Limpeza, geocodificação e filtragem de sequências
├── codon_analysis.py        # Análise de uso de códons e detecção de outliers
├── Nucleotide_analysis.py                    # Análise exploratória de caracteres inválidos
├── Align.py                 # Alinhamento das sequências por referência (paralelo)
│
├── ArqFasta/                # Sequências .fa baixadas (gerado pelo download.py)
├── FastaCombinados/         # FASTAs combinados por subtipo e segmento
├── Alinhamentos_Finais_RefBased/  # Resultados dos alinhamentos
│
├── cache_cidade_pais.json   # Cache de geocodificação
├── df_analise_N.pkl         # DataFrame intermediário salvo pelo preprocess.py
└── requirements.txt
```

---

## Dependências

```bash
pip install -r requirements.txt
```

**`requirements.txt`:**
```
biopython
pandas
numpy
matplotlib
seaborn
selenium
webdriver-manager
geopy
pycountry
pycountry-convert
googletrans==4.0.0-rc1
```

---

## Descrição dos Scripts

### `download.py`
Automatiza o download das sequências diretamente do [NCBI Influenza Virus Database](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi) usando Selenium. Itera sobre todos os subtipos H1–H18 × N1–N11 para os hospedeiros aviário, humano e suíno, baixando separadamente os segmentos HA e NA. Os arquivos são automaticamente renomeados no padrão `HxNx(hospedeiro)_SEGMENTO_.fa`.

### `preprocess.py`
Realiza o pipeline completo de limpeza e filtragem:
- Parsing dos identificadores FASTA com regex
- Correção de datas e remoção de registros inválidos
- Remoção de duplicatas
- Filtragem de stop codons incorretos (mantém apenas sequências iniciadas em ATG)
- **Geocodificação** das localidades para países e continentes (com cache JSON)
- Geração de gráficos cumulativos por data, país, continente e hospedeiro
- Exportação de FASTAs combinados por subtipo e segmento

### `Nucleotide_analysis.py`
Análise exploratória interativa de caracteres considerados inválidos (não-ACGT) nas sequências. Gera gráficos mostrando quantas sequências seriam removidas em diferentes limites percentuais de tolerância, auxiliando na escolha do threshold de filtragem.

### `codon_analysis.py`
Análise de uso de códons:
- Contagem e normalização por frequência relativa
- Boxplots com marcações de média e desvio padrão por segmento (HA e NA)
- Detecção de outliers via IQR por segmento
- Gráficos de correlação HA vs NA
- Comparação de médias de códons entre sequências com e sem outliers
- Exportação de tabelas CSV e FASTAs filtrados

### `Align.py`
Alinhamento par a par baseado em referência usando Biopython:
- Seleciona automaticamente a RefSeq específica para cada subtipo/segmento
- Usa RefSeq de fallback (H7N9) quando não há referência específica
- Processamento paralelizado

---

## Outputs Gerados

| Arquivo | Descrição |
|---------|-----------|
| `csvcerto2.csv` | DataFrame bruto após parsing |
| `df_analise_N.pkl` | DataFrame intermediário serializado |
| `codon_usage_table.csv` | Frequência relativa de cada códon por sequência |
| `estatisticas_codons_HA.csv` | Total, média e DP dos códons para HA |
| `estatisticas_codons_NA.csv` | Total, média e DP dos códons para NA |
| `Tabela_De_Outliers.csv` | Sequências com outliers e seus códons discrepantes |
| `df_com_out.csv` | IDs das sequências com outliers |
| `df_sem_out.csv` | IDs das sequências sem outliers |
| `df_sequencias_corretas.fa` | FASTA com sequências válidas após filtragem |
| `filtrado.fa` | FASTA com sequências identificadas como outliers |
| `FastaCombinados/` | FASTAs agrupados por subtipo e segmento |
| `Alinhamentos_Finais_RefBased/` | FASTAs alinhados por referência |

---

