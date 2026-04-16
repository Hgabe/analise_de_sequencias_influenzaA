from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import re
from multiprocessing import Pool, cpu_count # NOVO


PASTA_REFSEQS = "D:/sequencias para alinhamento" 

PASTA_ARQUIVOS_DE_DADOS = "D:/codigopython/progressdocodigo/FastaCombinados" 

REFSEQ_FALLBACK_HA = os.path.join(PASTA_REFSEQS, "(HA)(H7N9)NC_026425.fasta")
REFSEQ_FALLBACK_NA = os.path.join(PASTA_REFSEQS, "(NA)(H7N9)NC_026429.fasta")

# Mapeamento para escolher o fallback correto
MAPA_FALLBACK = {
    "HA": REFSEQ_FALLBACK_HA,
    "NA": REFSEQ_FALLBACK_NA
}

PASTA_SAIDA = "Alinhamentos_Finais_RefBased"

# ----------------------------#
#  PARÂMETROS DE ALINHAMENTO  #
# ----------------------------#

MATCH_SCORE = 4 #2 #teste2: 4
MISMATCH_SCORE = -4 
GAP_OPEN_PENALTY = -40  #teste2: -20 (dobrei de novo esse e o extend)
GAP_EXTEND_PENALTY = -12 #teste2: -6

# --------------------------#
#           FUNÇÕES         #
# --------------------------#

def extrair_subtipo_e_segmento(filename):
    nome_base = os.path.basename(filename)
    # REGEX : [Hh]\d{1,2} e [Nn]\d{1,2} para suportar H1, H10, H11, etc.
    match = re.search(r"([Hh]\d{1,2}[Nn]\d{1,2})\s*\((HA|NA)\)", nome_base, re.IGNORECASE)
    
    if match:
        subtipo = match.group(1).upper() 
        segmento = match.group(2).upper() 
        return segmento, subtipo
    return None, None


def encontrar_refseq_especifica(segmento, subtipo, pasta_refseqs):
    padrao = re.compile(f"\\({segmento}\\)\\({subtipo}\\)NC_.*\\.fasta", re.IGNORECASE)
    
    for filename in os.listdir(pasta_refseqs):
        if padrao.match(filename):
            return os.path.join(pasta_refseqs, filename)
    return None

def alinha_par_com_referencia(ref_seq, target_seq):
    alignments = pairwise2.align.globalms(
        ref_seq, 
        target_seq, 
        MATCH_SCORE, 
        MISMATCH_SCORE, 
        GAP_OPEN_PENALTY, 
        GAP_EXTEND_PENALTY
    )
    return alignments[0][1] 

def processar_arquivo_de_dados(arquivo_dados):
    nome_arquivo = os.path.basename(arquivo_dados)
    segmento, subtipo = extrair_subtipo_e_segmento(nome_arquivo) 
    
    if not segmento or not subtipo:
        return f"\nAVISO: Padrão não reconhecido no arquivo {nome_arquivo}. Pulando."

    print(f"\n--- Processando {subtipo} ({segmento}) ---")
    
    try:
        ref_path = encontrar_refseq_especifica(segmento, subtipo, PASTA_REFSEQS)
        
        if ref_path:
            print(f"  -> Usando RefSeq Específica: {os.path.basename(ref_path)}")
        else:
            if segmento in MAPA_FALLBACK:
                ref_path = MAPA_FALLBACK[segmento]
                print(f"  -> RefSeq Específica ({subtipo}) não encontrada. Usando Fallback H7N9 ({segmento}).")
            else:
                return f"ERRO FATAL: Não há RefSeq de fallback configurada para o segmento '{segmento}'."

        if not os.path.exists(ref_path):
             return f"ERRO FATAL: Arquivo de referência não encontrado em {ref_path}"
        
   
        ref_record = SeqIO.read(ref_path, "fasta")
        ref_seq = str(ref_record.seq)
        arquivo_saida = os.path.join(PASTA_SAIDA, f"{subtipo}_{segmento}_ALINHADO_RefBased.fasta")

        SeqIO.write(ref_record, arquivo_saida, "fasta")
        sequencias_alinhadas = []

        for record in SeqIO.parse(arquivo_dados, "fasta"):
            if record.id != ref_record.id: 
                aligned_seq = alinha_par_com_referencia(ref_seq, str(record.seq))
                
                aligned_record = SeqRecord(
                    Seq(aligned_seq), 
                    id=record.id, 
                    name=record.name, 
                    description=record.description
                )
                sequencias_alinhadas.append(aligned_record)
        
        with open(arquivo_saida, "a") as output_handle:
            SeqIO.write(sequencias_alinhadas, output_handle, "fasta")

        return f"Alinhamento para {subtipo} ({segmento}) salvo em {arquivo_saida}"
        
    except Exception as e:
        return f"ERRO no processamento de {arquivo_dados}: {e}"

# ----------------------------------#
# PROCESSO PRINCIPAL (PARALELIZADO) #
# ----------------------------------#

if __name__ == '__main__':
    os.makedirs(PASTA_SAIDA, exist_ok=True)
    arquivos_de_dados = [
        os.path.join(PASTA_ARQUIVOS_DE_DADOS, f) 
        for f in os.listdir(PASTA_ARQUIVOS_DE_DADOS) 
        if f.endswith(('.fasta', '.fa'))
    ]
    NUM_PROCESSOS = max(1, cpu_count() // 2) 
    print(f"\nUsando {NUM_PROCESSOS} processos paralelos")
    
    with Pool(NUM_PROCESSOS) as pool:
        resultados = pool.map(processar_arquivo_de_dados, arquivos_de_dados)
    
    for resultado in resultados:
        print(resultado)

    print("\nProcesso de Alinhamento Baseado em Referência Concluído")