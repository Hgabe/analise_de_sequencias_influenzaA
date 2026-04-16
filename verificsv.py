import os
import re
import json
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter

import pycountry
import pycountry_convert as pc
from geopy.geocoders import Nominatim
from geopy.exc import GeocoderTimedOut, GeocoderServiceError
from googletrans import Translator

# ---------------------------------------------------------------------------- #
#                              CONFIGURAÇÕES GLOBAIS                           #
# ---------------------------------------------------------------------------- #

HERE          = os.path.abspath(os.path.dirname(__file__))
DOWNLOAD_PATH = os.path.join(HERE, 'ArqFasta')
COM_CODONS    = True
CACHE_FILE    = "cache_cidade_pais.json"
STOP_CODONS   = ["TAA", "TGA", "TAG"]

IDENTIFIER_PATTERN = re.compile(
    r"(?P<ID>[^\s]+)\s(?P<Type>[A-Z])\/(?P<Species>[^/]+(?:\s[^/]+)*)\/(?P<Location>[^\/]+)\/(?P<Data>[^/]+)\/(?P<Sequence_type>[^\s]+)\s"
)

# Mapeamento manual de cidades/regiões para países
MOCK_LLM_RESPONSES = {
    "Hunan": "China", "Tsukuba": "Japan", "Korea": "Korea, Republic of", "Hokkaido": "Japan", "Zhejiang": "China",
    "Jiangxi": "China", "Northern Shoveler": None, "Fujian": "China", "Saga": "Japan", "Hubei": "China",
    "Kumamoto": "Japan", "Shimane": "Japan", "Shandong": "China", "Okayama": "Japan", "Shanxi": "China",
    "Guizhou": "China", "Beijing": "China", "Hebei": "China", "Guangdong": "China", "Shanghai": "China",
    "Jiangsu": "China", "Miyagi": "Japan", "Nanjing": "China", "Yangzhou": "China", "Xuyi": "China",
    "Hongze": "China", "Zhalong": "China", "Chiba": "Japan", "Yunnan": "China", "Sanjiang": "China",
    "Caizi Lake": "China", "Shengjin Lake": "China", "Wuhan": "China", "Dongting": "China", "Weihai": "China",
    "CHBZ": None, "Liaoning": "China", "Dalian": "China", "Novomychalivka": "Ukraine", "Gurjev": "Russia",
    "Anhui": "China", "Aomori": "Japan", "Nanchang": "China", "Yasnopolyanske": "Ukraine", "Aichi": "Japan",
    "Ahrenswolde": "Germany", "Egersheim": "Germany", "Eire": "Ireland", "Finnistere": "France",
    "Grossharie": "Germany", "Guangxi": "China", "Gueterlsoh": "Germany", "Heilongjiang": "China",
    "Henan": "China", "Hiroshima": "Japan", "Ibaraki": "Japan", "Jilin": "China", "Kagoshima": "Japan",
    "Krerensen": "Germany", "Kyoto": "Japan", "Laiwu": "China", "Minnoosta": "United States",
    "Narita": "Japan", "Niigata": "Japan", "Ningjin": "China", "Okinawa": "Japan", "Ortensburg": "Germany",
    "Osaka": "Japan", "Osaka-c": "Japan", "Pingtung": "Taiwan", "Raesefeld": "Germany", "Shaanxi": "China",
    "Sichuan": "China", "SouthDakota": "United States", "an": None, "Taichung": "Taiwan", "Tianjin": "China",
    "Tottori": "Japan", "Wetedorf": "Germany", "Yamagata": "Japan", "Zhucheng": "China", "Akita": "Japan",
    "Albany": "United States", "Amagasaki": "Japan", "Azarbayejan_Sharghi": "Iran", "Baofeng": "China",
    "Bayanulgii": "Mongolia", "Blagovechensk": "Russia", "CHR": None, "Cameron": "United States",
    "Changchun": "China", "Changsha": "China", "Chengdu": "China", "Dongcheng": "China",
    "EastBaltimore": "United States", "FuZhou": "China", "Fukushima": "Japan", "Fuzhou": "China",
    "Guangzhou": "China", "GuangzhouSB": "China", "Gunma": "Japan", "Gyeonbuk": "Korea, Republic of",
    "Haishu": "China", "Hangzhou": "China", "Himeji": "Japan", "HuZhou": "China", "Hualong": "China",
    "Hunan Changsha": "China", "Hunan Furong": "China", "Hunan Kaifu": "China", "Hyogo": "Japan",
    "Inner Mongolia": "Mongolia", "Iwate": "Japan", "Jiangsubeitang": "China", "Jiangsudanyang": "China",
    "Jiangsugaoyou": "China", "Jiangsugusu": "China", "Jiangsuhailing": "China", "Jiangsunanchang": "China",
    "Jiangsunanjing": "China", "Jiangsupizhou": "China", "Jiangsuqinhuai": "China", "Jiangsuquanshan": "China",
    "Jiangsusucheng": "China", "Jiangsutinghu": "China", "Jiangyin": "China", "Jinhua": "China",
    "Jinshui": "China", "Kanagawa": "Japan", "Kawasaki": "Japan", "Keelung": "Taiwan",
    "Kobe": "Japan", "Kowloon": "China", "Linkou": "Taiwan", "Lishui": "China", "LongYan": "China",
    "Mie": "Japan", "Morioka": "Japan", "Moscowv": "Russia", "Nagano": "Japan", "Nagasaki": "Japan",
    "NanChang": "China", "NanPing": "China", "Neimenggu": "China", "NewJersey": "United States",
    "NingDe": "China", "Para": "Brazil", "PuTian": "China", "Qingdao": "China", "Qingfeng": "China",
    "QuanZhou": "China", "Ruyang": "China", "Ruzhou": "China", "Saitama": "Japan", "Sakai": "Japan",
    "SanMing": "China", "Sapporo": "Japan", "Sendai": "Japan", "Shantou": "China", "Shengzhen": "China",
    "Shenzhen": "China", "Shiga": "Japan", "Shizuoka": "Japan", "StEtienne": "France", "Suita": "Japan",
    "Taipei": "Taiwan", "Taizhou": "China", "TayNguyen": "Vietnam", "TayNinh": "Vietnam", "ThaiBinh": "Vietnam",
    "Tianjinhedong": "China", "Tianjinheping": "China", "Tianjinhongqiao": "China", "Tianjinjinnan": "China",
    "Tianjinnankai": "China", "Tianjininghe": "China", "Tianjintanggu": "China", "Tientsin": "China",
    "Tochigi": "Japan", "Tokushima": "Japan", "Tokyo": "Japan", "Tumen": "China", "Unmugobi": "Mongolia",
    "Utsunomiya": "Japan", "Viet nam": "Vietnam", "Vladivistok": "Russia", "Wakayama": "Japan",
    "Weidong": "China", "Weidu": "China", "XiaMen": "China", "Xiamen": "China", "Xian": "China",
    "Xigong": "China", "Xuchang": "China", "Yamaguchi": "Japan", "Yokohama": "Japan", "Yokosuka": "Japan",
    "Ytml": None, "Zanjann": "Iran", "ZhangZhou": "China", "Zhongyuan": "China", "Zhoushan": "China",
    "BaRiaVungTau": "Vietnam", "BinhDuong": "Vietnam", "Ehime": "Japan", "NamDinh": "Vietnam",
    "TienGiang": "Vietnam", "Zhaotong": "China", "Yuhuan": "China", "Hainan": "China", "Yunan": "China",
    "Binh Doung": "Vietnam", "Gescheu": "Germany", "Goldenbeck": "Germany", "Langefoerden": "Germany",
    "Miyazaki": "Japan", "Salingberg": "Germany", "Spelk": "Germany", "Tainan": "Taiwan", "spain": "Spain",
    "Wenzhou": "China", "Praimoric": "Russia", "Wuxi": "China", "Adachi": "Japan", "Izumi": "Japan",
    "Zhang": "China", "Longquan": "China", "Xianghai": "China", "Wuerusreuth": "Germany", "Foshan": "China",
    "Ganzhou": "China", "GuangXi": "China", "Jiangshu": "China", "turkey": "Turkey", "Boitzenborstel": "Germany",
    "Harmstrub": "Germany", "HuNan": "China", "Karrenzien": "Germany", "Missourri": "United States",
    "Obihiro": "Japan", "BacGiang": "Vietnam", "Beijingdongcheng": "China", "Beijingfengtai": "China",
    "Beijinghuairou": "China", "Beijingxicheng": "China", "Columbia": "United States", "Fukuoka": "Japan",
    "GREECE": "Greece", "Gifu": "Japan", "Guandong": "China", "HONG KONG": "China", "Harbin": "China",
    "HoaBinh": "Vietnam", "Kong SAR": "China", "Illiniois": "United States", "Ishikawa": "Japan",
    "Kagawa": "Japan", "PERU": "Peru", "PHILIPPINES": "Philippines", "Qinghai": "China", "SINGAPORE": "Singapore",
    "SWITZERLAND": "Switzerland", "Suzhou": "China", "Voroneg": "Russia", "Wujiaqu": "China",
    "Zhongshan": "China", "Eastern China": "China", "Zhuanghe": "China", "winged Teal": None, "Altai": "Russia",
    "SanJiang": "China", "Ebinur Lake": "China", "Dongguan": "China", "Huizhou": "China", " Mongolia": "Mongolia",
    "Nunavet": "Canada", "Chibi": "China", "HuBei": "China", "Khurkh river": "Mongolia", "PoyangLake": "China",
    "Poyang Lake": "China", "ZhaLong": "China", "Chongqing": "China", "Deli Derdang": "Indonesia",
    "Murao Jambi": "Indonesia", "Papua": "Indonesia", "Anyang": "China", "Huadong": "China", "Taoyuan": "Taiwan",
    "Tibet": "China", "AnNing": "China", "CentralJava": "Indonesia", "DaLi": "China", "Danyang": "China",
    "Dongtai": "China", "EastJava": "Indonesia", "EastKalimantan": "Indonesia", "Gansu": "China", "Gaoyou": "China",
    "Huabei": "China", "Kohn Kaen": "Thailand", "Lao": "Laos", "Menofia": "Egypt", "Nakornsawan": "Thailand",
    "Nara": "Japan", "Ningxia": "China", "Oita": "Japan", "Qalubia": "Egypt", "Sharkeya": "Egypt",
    "Sutiakhali": "Bangladesh", "TongHai": "China", "Tonghai": "China", "WestJava": "Indonesia",
    "Xinjiang": "China", "Xuzhou": "China", "Zhaozhuang": "China", "VinhLong": "Vietnam",
    "'Ivoire": "Côte d'Ivoire", "Phuyen": "Vietnam", "Quangngai": "Vietnam", "Thai Ninh": "Vietnam",
    "Vietnam BacLiu": "Vietnam", "Vietnam NamDinh": "Vietnam", "Krasnoozerskoye": "Russia",
    "headed gull": None, "Panjing": "China", "Long An": "Vietnam", "Muaraenim": "Indonesia", "Toyama": "Japan",
    "Chakshahzad": "Pakistan", "Phathumthani": "Thailand", "Yogjakarta": "Indonesia", "Cygnus olor": None,
    "Aguascallientes": "Mexico", "Changhua": "Taiwan", "Chiayi": "Taiwan", "Chiping": "China", "Hsinchu": "Taiwan",
    "Kaohsiung": "Taiwan", "Miaoli": "Taiwan", "Taipei City": "Taiwan", "Queretaroa": "Mexico",
    "Shijiazhuang": "China", "Yunlin": "Taiwan", "Jiang Xi": "China", "Nantou": "Taiwan", "Yilan": "Taiwan",
    "Djankoy": "Ukraine", "bird feces": None, "Taishun": "China", "Khunt lake": "Mongolia",
    "Sanmenxia": "China", "kumamoto": "Japan", "Egypt ": "Egypt", "Yuhang": "China", "Nanji": "China",
    "San Jiang": "China", "Hei Longjiang": "China", "Tumuji": "China", "Antarctic": "Antarctica",
    "XiangHai": "China", "Hongkong": "China", "Maduo": "China", "Gonghe": "China", "Gangcha": "China",
    "Hadoree": "Korea, Republic of", "winged teal": None, "GSC_chicken_B": None, "Huzhou": "China",
    "Shaoxing": "China", "Dalai lake": "Mongolia", "Fukui": "Japan", "Nizhnegorsk": "Ukraine",
    "Novomyhalivka": "Ukraine", "Heinan": "China", "Inner_Mongolia": "China", "Jiaxing": "China",
    "Quzhou": "China", "Rizhao": "China", "Tianjing": "China", "Zhangzhou": "China", "Sunan": "China",
    "eastern China": "China", "Kunming": "China", "Minhang": "China", "Qingyuan": "China", "WuXi": "China",
    "Zhenjiang": "China", "wuhu": "China", "Baicheng": "China", "Gangxi": "China", "Daye": "China",
    "EGYPT": "Egypt", "FuJian": "China", "HeBei": "China", "HeNan": "China", "Heibei": "China",
    "Heiilongjiang": "China", "HongKong": "China", "Huaian": "China", "Jiande": "China", "JinShui": "China",
    "Jinan": "China", "Jingmen": "China", "Juxian": "China", "Kasoor": "Pakistan", "Mudanjiang": "China",
    "Nantong": "China", "Ningbo": "China", "Sandong": "China", "ShanXi": "China", "ShangDong": "China",
    "Shangdong": "China", "Shaoguan": "China", "Shenyang": "China", "Suqian": "China", "Wujin": "China",
    "XinjiangBaicheng": "China", "Xuancheng": "China", "Yancheng": "China", "Yantai": "China",
    "ShanDong": "China", "InnerMongolia": "China", "Lengshuitan": "China", "Changzhou": "China",
    "Azarbayejan-Sharghi": "Iran",
}

CONTINENTE_MAP = {
    'NA': 'América do Norte', 'SA': 'América do Sul', 'AS': 'Ásia',
    'EU': 'Europa', 'AF': 'África', 'OC': 'Oceania', 'AN': 'Antártica'
}


# ---------------------------------------------------------------------------- #
#                              LEITURA DOS DADOS                               #
# ---------------------------------------------------------------------------- #

def carregar_dados(download_path: str) -> list:
    data = []
    for file_name in os.listdir(download_path):
        if not file_name.endswith(".fa"):
            continue

        file_path = os.path.join(download_path, file_name)

        class_name = re.search(r"\((.*?)\)", file_name)
        class_name = class_name.group(1) if class_name else "Unknown"

        segment = re.search(r"_(.*?)_", file_name)
        segment = segment.group(1) if segment else "Unknown"

        subtype = re.search(r"^(H\d+N\d+)", file_name)
        subtype = subtype.group(1) if subtype else "Unknown"

        with open(file_path, "r") as fasta_file:
            identifier, sequence = "", ""
            for line in fasta_file:
                line = line.strip()
                if line.startswith(">"):
                    if identifier and sequence:
                        data.append([identifier, sequence, class_name, segment, subtype])
                    identifier = line[1:]
                    sequence = ""
                else:
                    sequence += line
            if identifier and sequence:
                data.append([identifier, sequence, class_name, segment, subtype])

    return data


# ---------------------------------------------------------------------------- #
#                         CONSTRUÇÃO E LIMPEZA DO DATAFRAME                    #
# ---------------------------------------------------------------------------- #

def construir_dataframe(data: list) -> pd.DataFrame:
    data_cleaned = []
    for entry in data:
        match = IDENTIFIER_PATTERN.match(entry[0])
        if match:
            data_cleaned.append([
                match.group("ID"), match.group("Type"), match.group("Species"),
                match.group("Location"), match.group("Data"), match.group("Sequence_type"),
                entry[1], entry[2], entry[3], entry[4]
            ])

    return pd.DataFrame(data_cleaned, columns=[
        "ID", "Type", "Species", "Location", "Sequence_type", "Data",
        "Sequence", "Class", "Segment", "Subtype"
    ])


def corrigir_colunas(df: pd.DataFrame) -> pd.DataFrame:
    df[["Data", "Sequence_type", "Location", "Species"]] = df.apply(
        lambda row: pd.Series([
            row["Sequence_type"] if "/" in row["Data"] else row["Data"],
            row["Location"]      if "/" in row["Data"] else row["Sequence_type"],
            row["Species"]       if "/" in row["Data"] else row["Location"],
            np.nan               if "/" in row["Data"] else row["Species"]
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


def remover_duplicatas(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(subset=["ID", "Sequence"], keep='last')


def filtrar_location_invalida(df: pd.DataFrame) -> pd.DataFrame:
    return df[~(df["Location"].str.contains(r"\d", regex=True) | (df["Location"].str.len() <= 2))]


def limpar_location(df: pd.DataFrame) -> pd.DataFrame:
    df['Location'] = df['Location'].str.replace(r'-\s?[Cc]$',  '', regex=True, flags=re.IGNORECASE)
    df['Location'] = df['Location'].str.replace(r'-\s?SB$',    '', regex=True, flags=re.IGNORECASE)
    df['Location'] = df['Location'].str.strip()
    return df


def filtrar_sequencias_invalidas(df: pd.DataFrame, threshold: float = 0.01) -> pd.DataFrame:
    def pct_invalidos(seq):
        total = len(seq)
        invalidos = sum(1 for c in seq if c not in "AGCT")
        return invalidos / total if total > 0 else 1

    df["Invalid_pct"] = df["Sequence"].apply(pct_invalidos)
    df = df[df["Invalid_pct"] < threshold]
    return df


def filtrar_stop_codons(df: pd.DataFrame) -> pd.DataFrame:
    def count_stop_codons(seq):
        return sum(1 for i in range(0, len(seq) - 2, 3) if seq[i:i+3] in STOP_CODONS)

    df["StopCodons_count"] = df["Sequence"].apply(count_stop_codons)
    mask_atg       = df["Sequence"].str.startswith("ATG")
    mask_sem_stop  = df["StopCodons_count"] == 0
    mask_com_stop  = (df["StopCodons_count"] == 1) & df["Sequence"].str.endswith(("TAA", "TAG", "TGA"))
    df = df[mask_atg & (mask_sem_stop | mask_com_stop)].copy()
    df.drop(columns=["StopCodons_count"], inplace=True)
    return df


# ---------------------------------------------------------------------------- #
#                           GEOCODIFICAÇÃO E PAÍSES                            #
# ---------------------------------------------------------------------------- #

def carregar_cache(cache_file: str) -> dict:
    if os.path.exists(cache_file) and os.path.getsize(cache_file) > 0:
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"Erro ao carregar cache, recriando: {e}")
    return {}


def salvar_cache(cache: dict, cache_file: str):
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)


def get_country_name(country_code: str) -> str | None:
    try:
        country = pycountry.countries.get(alpha_2=country_code)
        return country.name if country else None
    except Exception:
        return None


def traduzir_nome(nome: str, cache: dict, translator) -> str:
    if nome in cache:
        return cache[nome]
    try:
        translated = translator.translate(nome, dest='en').text
        cache[nome] = translated
        print(f"Traduzido '{nome}' -> '{translated}'")
        return translated
    except Exception:
        return nome


def geocode_nominatim(cidade: str, geolocator, translator, traducao_cache: dict, use_translation: bool = False, timeout: int = 10) -> str | None:
    if use_translation:
        cidade = traduzir_nome(cidade, traducao_cache, translator)
    try:
        loc = geolocator.geocode(cidade, exactly_one=True, timeout=timeout, addressdetails=True)
        if loc and loc.raw and 'address' in loc.raw:
            address = loc.raw['address']
            country_code = address.get('country_code')
            if country_code:
                return get_country_name(country_code.upper())
            country_name = address.get('country')
            if country_name:
                try:
                    code = pc.country_name_to_country_alpha2(country_name, cn_name_format="default")
                    return get_country_name(code)
                except Exception:
                    return country_name
    except Exception:
        return None
    return None


def cidade_para_pais_master(cidade: str, geolocator, translator, traducao_cache: dict, timeout: int = 10) -> str | None:
    pais = geocode_nominatim(cidade, geolocator, translator, traducao_cache, use_translation=False, timeout=timeout)
    if pais:
        print(f"Nominatim direto: '{cidade}' -> '{pais}'")
        return pais

    pais = geocode_nominatim(cidade, geolocator, translator, traducao_cache, use_translation=True, timeout=timeout)
    if pais:
        print(f"Nominatim c/ tradução: '{cidade}' -> '{pais}'")
        return pais

    pais = MOCK_LLM_RESPONSES.get(cidade)
    if pais:
        print(f"Dicionário manual: '{cidade}' -> '{pais}'")
    else:
        print(f"Nenhum resultado para '{cidade}'")
    return pais


def geocodificar_cidades(df: pd.DataFrame, cache_file: str, save_every: int = 50, max_tries: int = 3) -> pd.DataFrame:
    geolocator   = Nominatim(user_agent="city-to-country-script-v5-final")
    translator   = Translator()
    traducao_cache = {}
    cidade_para_pais_dict = carregar_cache(cache_file)

    cidades_unicas = df["Location"].dropna().unique()
    to_check = [c for c in cidades_unicas
                if c not in cidade_para_pais_dict or cidade_para_pais_dict[c] is None]

    for i, cidade in enumerate(to_check, 1):
        cidade_str = str(cidade)
        if cidade_str.isdigit() or len(cidade_str) <= 2:
            cidade_para_pais_dict[cidade] = None
            continue

        print(f"[{i}/{len(to_check)}] Consultando: {cidade_str}")
        success = False
        for attempt in range(1, max_tries + 1):
            try:
                pais = cidade_para_pais_master(cidade_str, geolocator, translator, traducao_cache)
                cidade_para_pais_dict[cidade] = pais
                success = True
                break
            except GeocoderTimedOut:
                wait = 2 ** attempt
                print(f"Timeout (tentativa {attempt}). Backoff {wait}s...")
                time.sleep(wait)
            except GeocoderServiceError as e:
                print(f"Erro de serviço: {e}. Aguardando 10s...")
                time.sleep(10)
            except Exception as e:
                print(f"Erro inesperado em '{cidade_str}': {e}")
                break

        if not success:
            cidade_para_pais_dict[cidade] = None

        time.sleep(1)

        if i % save_every == 0:
            salvar_cache(cidade_para_pais_dict, cache_file)
            print(f"Progresso salvo ({i} consultas).")

    salvar_cache(cidade_para_pais_dict, cache_file)
    print("Cache atualizado e salvo.")

    df["Country"] = df["Location"].map(cidade_para_pais_dict)
    return df


def adicionar_continente(df: pd.DataFrame) -> pd.DataFrame:
    continente_cache = {}

    def pais_para_continente(pais):
        if pais in continente_cache:
            return continente_cache[pais]
        try:
            code = pc.country_name_to_country_alpha2(pais, cn_name_format="default")
            cont = pc.country_alpha2_to_continent_code(code)
            resultado = CONTINENTE_MAP.get(cont, 'Outros')
        except Exception:
            resultado = 'Outros'
        continente_cache[pais] = resultado
        return resultado

    df["Continente"] = df["Country"].apply(pais_para_continente)
    return df


# ---------------------------------------------------------------------------- #
#                               ANÁLISE DE CÓDONS                              #
# ---------------------------------------------------------------------------- #

def contar_codons(seq: str) -> dict:
    return dict(Counter(seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3))


def gerar_tabela_codons(df: pd.DataFrame, normalizar: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    codon_dicts = df["Sequence"].apply(contar_codons)
    codon_df = pd.DataFrame.from_records(codon_dicts).fillna(0).astype(int).reset_index(drop=True)

    codon_df_norm = codon_df.copy().drop(columns=STOP_CODONS, errors='ignore')
    if normalizar:
        codons_count = df["Codons_count"].reset_index(drop=True)
        codon_df_norm = codon_df_norm.div(codons_count, axis=0)

    return codon_df, codon_df_norm


def relatorio_stop_codons(df1: pd.DataFrame, codon_df_norm: pd.DataFrame):
    def count_stop(seq):
        return sum(1 for i in range(0, len(seq) - 2, 3) if seq[i:i+3] in STOP_CODONS)

    df1_copy = df1.copy()
    df1_copy["StopCodons_count"] = df1_copy["Sequence"].apply(count_stop)

    zero   = len(df1_copy[df1_copy["StopCodons_count"] == 0])
    um     = len(df1_copy[df1_copy["StopCodons_count"] == 1])
    mais   = df1_copy[df1_copy["StopCodons_count"] > 1]["StopCodons_count"].value_counts().sort_index()
    total_mais = sum(mais)

    print(f"Sequências com 0 códons de parada:       {zero}")
    print(f"Sequências com exatamente 1 stop codon:  {um}")
    print(f"Sequências com mais de 1 stop codon:     {total_mais}")
    print("\nDetalhamento (mais de 1 stop codon):")
    print(mais.to_string())

    df_one = df1_copy[df1_copy["StopCodons_count"] == 1]
    print(f"\nTerminadas em TAA: {df_one['Sequence'].str.endswith('TAA').sum()}")
    print(f"Terminadas em TAG: {df_one['Sequence'].str.endswith('TAG').sum()}")
    print(f"Terminadas em TGA: {df_one['Sequence'].str.endswith('TGA').sum()}")


def detectar_outliers_por_segmento(df1: pd.DataFrame, codon_df_norm: pd.DataFrame) -> tuple:
    def iqr_mask(grupo):
        Q1, Q3 = grupo.quantile(0.25), grupo.quantile(0.75)
        IQR = Q3 - Q1
        return (grupo < (Q1 - 1.5 * IQR)) | (grupo > (Q3 + 1.5 * IQR))

    outlier_mask = codon_df_norm.groupby(df1["Segment"]).transform(iqr_mask)
    outlier_df   = codon_df_norm[outlier_mask]

    tabela_raw = pd.concat([df1["ID"].reset_index(drop=True), outlier_df.reset_index(drop=True)], axis=1)

    def formatar_linha(row):
        return [f"{val}({col})" for col, val in row.items() if pd.notna(val)]

    linhas = tabela_raw.apply(formatar_linha, axis=1)
    Tabela_De_Outliers = pd.DataFrame({"Valores": [", ".join(l) for l in linhas]})
    Tabela_De_Outliers = Tabela_De_Outliers[Tabela_De_Outliers["Valores"].str.count(r"\(") >= 2]

    df_out_temp = pd.concat([df1["ID"].reset_index(drop=True), outlier_df.reset_index(drop=True)], axis=1)
    df_out_temp["tem_out"] = (~df_out_temp.drop(columns="ID").isna().all(axis=1)).astype(int)
    df_com_out = df_out_temp[df_out_temp["tem_out"] == 1].drop(columns="tem_out")
    df_sem_out = df_out_temp[df_out_temp["tem_out"] == 0].drop(columns="tem_out")

    return Tabela_De_Outliers, df_com_out, df_sem_out, outlier_mask


def codons_por_grupo(df1: pd.DataFrame, codon_df: pd.DataFrame, ids: pd.Series, segment: str) -> pd.DataFrame:
    df_filtrado = df1[df1["ID"].isin(ids)].drop_duplicates(subset=["ID"]).reset_index(drop=True)
    resultado   = pd.concat([df_filtrado[["ID", "Segment"]], codon_df.reset_index(drop=True)], axis=1)
    resultado   = resultado.dropna(subset=["ID"])
    return resultado[resultado["Segment"] == segment]


# ---------------------------------------------------------------------------- #
#                               VISUALIZAÇÕES                                  #
# ---------------------------------------------------------------------------- #

def plotar_countplot(df: pd.DataFrame, x: str, titulo: str, xlabel: str, ylabel: str, log: bool = True, figsize: tuple = (8, 5)):
    plt.figure(figsize=figsize)
    sns.countplot(x=x, data=df, color='gray', order=df[x].value_counts().index)
    plt.title(titulo, fontsize=16)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    if log:
        plt.yscale("log")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def plotar_countplot_por_segmento(df: pd.DataFrame, segmentos: list):
    for seg in segmentos:
        df_seg = df[df['Segment'] == seg]
        contagens = df_seg['Class'].value_counts()
        print(f'\nContagem de Sequências por Hospedeiro para o gene {seg}')
        print(contagens)
        plt.figure(figsize=(8, 5))
        sns.countplot(x='Class', color='gray', data=df_seg, order=contagens.index)
        plt.title(f'Contagem de Sequências por Hospedeiro — Gene {seg}', fontsize=14)
        plt.xlabel('Hospedeiro', fontsize=12)
        plt.ylabel('Número de Sequências', fontsize=12)
        plt.yscale("log")
        plt.xticks(rotation=90)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.show()


def plotar_cumulativo(df: pd.DataFrame, segmento: str, groupby_col: str, col_label: str, titulo: str):
    df_filtrado = df[df["Segment"] == segmento]
    df_counts   = df_filtrado.groupby(["Data", groupby_col]).size().unstack(fill_value=0)
    df_cumsum   = df_counts.cumsum().apply(pd.to_numeric, errors="coerce").fillna(0)

    df_cumsum.plot(kind="line", stacked=False, title=titulo)
    plt.xticks(ticks=df_cumsum.index, labels=df_cumsum.index,
               rotation=90, ha='center', fontsize=7)
    plt.yscale("log")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel("Data de Coleta")
    plt.ylabel("Quantidade Acumulada (log)")
    plt.tight_layout()
    plt.show()


def plotar_cumulativo_geral(df: pd.DataFrame):
    df_counts = df.groupby(["Data", "Segment"]).size().unstack(fill_value=0)
    df_cumsum = df_counts.cumsum()
    df_cumsum.plot(kind="line", stacked=False, color=["blue", "red"],
                   title="Curva Cumulativa de Sequências por Data de Coleta e Segmento")
    plt.xticks(ticks=df_cumsum.index, labels=df_cumsum.index,
               rotation=90, ha='center', fontsize=7)
    plt.yscale("log")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel("Data de Coleta")
    plt.ylabel("Quantidade Acumulada (log)")
    plt.tight_layout()
    plt.show()


def plotar_boxplot_codons(codon_dfs: pd.DataFrame, codon_df: pd.DataFrame, segment: str, download_path: str):
    estatisticas = pd.DataFrame({
        "total":        codon_df.sum(),
        "Media":        codon_df.mean(),
        "DesvioPadrao": codon_df.std()
    })
    estatisticas.to_csv(os.path.join(download_path, f"estatisticas_codons_{segment}.csv"))

    codon_df_sorted = codon_dfs[sorted(codon_df.columns)]
    mean_values = codon_df_sorted.mean()
    std_values  = codon_df_sorted.std()

    sns.boxplot(data=codon_df_sorted, color="lightgray",
                boxprops=dict(facecolor='gray'), showfliers=False)

    for i, (mean, std) in enumerate(zip(mean_values, std_values)):
        plt.text(i, std,  '__', ha='center', va='center', color='blue', fontsize=10)
        plt.text(i, mean, '__', ha='center', va='center', color='red',  fontsize=10)

    media_patch   = mpatches.Patch(color='red',   label='Média')
    dp_patch      = mpatches.Patch(color='blue',  label='Desvio Padrão')
    mediana_patch = mpatches.Patch(color='black', label='Mediana')
    plt.legend(handles=[media_patch, dp_patch, mediana_patch],
               loc='upper right', bbox_to_anchor=(1, 1))

    plt.title(f"Box Plot com Média e Desvio Padrão dos Códons — Influenza A ({segment})")
    plt.xlabel("Códons")
    plt.ylabel("Frequência Relativa")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def plotar_boxplot_simples(codon_df: pd.DataFrame, titulo: str):
    codon_df_sorted = codon_df[sorted(codon_df.columns)]
    sns.boxplot(data=codon_df_sorted, color="lightgray",
                boxprops=dict(facecolor='gray'), showfliers=False)
    plt.title(titulo)
    plt.xlabel("Códons")
    plt.ylabel("Contagem")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def plotar_comparacao_barras(comparacao: pd.DataFrame, segment: str):
    comparacao.plot(kind="bar", color=["gray", "black"], alpha=0.8)
    plt.xlabel("Códons")
    plt.ylabel("Contagem Média")
    plt.title(f"Comparação de Média de Códons com e sem Outliers — Gene {segment}")
    plt.xticks(rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plotar_dispersao_medias(comparacao: pd.DataFrame, segment: str):
    codons  = comparacao.index
    x_pos   = range(len(codons))
    col_sem = "Não outlier"
    col_com = "Outlier"

    plt.figure(figsize=(12, 6))
    plt.scatter(x_pos, comparacao[col_sem], color="gray",  label=col_sem, marker="o", s=50)
    plt.scatter(x_pos, comparacao[col_com], color="black", label=col_com, marker="x", s=50)
    plt.xlabel("Códons")
    plt.ylabel("Contagem Média")
    plt.title(f"Dispersão de Médias com e sem Outliers — Gene {segment}")
    plt.xticks(x_pos, codons, rotation=90)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()


def plotar_correlacao_ha_na(df_sem: pd.DataFrame, df_com: pd.DataFrame):
    for grupo_label, df_ha, df_na in [
        ("Sem Outliers",
         df_sem[df_sem["Segment"] == "HA"].drop(columns=["ID", "Segment"]),
         df_sem[df_sem["Segment"] == "NA"].drop(columns=["ID", "Segment"])),
        ("Com Outliers",
         df_com[df_com["Segment"] == "HA"].drop(columns=["ID", "Segment"]),
         df_com[df_com["Segment"] == "NA"].drop(columns=["ID", "Segment"])),
    ]:
        med_ha = df_ha.median(numeric_only=True)
        med_na = df_na.median(numeric_only=True)
        df_med = pd.DataFrame({"HA": med_ha, "NA": med_na}).sort_values("HA")

        plt.figure(figsize=(12, 8))
        sns.scatterplot(x="HA", y="NA", data=df_med)
        plt.plot([0, df_med["HA"].max()], [0, df_med["HA"].max()],
                 'r--', label='y=x (Correlação Perfeita)')
        for i, txt in enumerate(df_med.index):
            plt.annotate(txt, (df_med["HA"].iloc[i], df_med["NA"].iloc[i]),
                         xytext=(-4, 4), textcoords='offset points', fontsize=8, alpha=0.7)
        plt.title(f"Correlação do Uso de Códons (HA vs NA) — {grupo_label}")
        plt.xlabel("Frequência Relativa Mediana (Gene HA)")
        plt.ylabel("Frequência Relativa Mediana (Gene NA)")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend()
        plt.show()


def plotar_impacto_outliers(med_sem: pd.Series, med_com: pd.Series, segment: str):
    df_comp = pd.DataFrame({
        f'{segment} Sem Outlier': med_sem,
        f'{segment} Com Outlier': med_com
    }).sort_values(f'{segment} Sem Outlier')

    col_sem = f'{segment} Sem Outlier'
    col_com = f'{segment} Com Outlier'
    max_val = max(df_comp[col_sem].max(), df_comp[col_com].max())

    plt.figure(figsize=(12, 8))
    sns.scatterplot(x=col_sem, y=col_com, data=df_comp)
    plt.plot([0, max_val], [0, max_val], 'r--', label='y=x (Nenhum Impacto do Outlier)')
    for i, txt in enumerate(df_comp.index):
        plt.annotate(txt, (df_comp[col_sem].iloc[i], df_comp[col_com].iloc[i]),
                     xytext=(-4, 4), textcoords='offset points', fontsize=8, alpha=0.7)
    plt.title(f"Impacto de Outliers na Mediana de Frequência de Códons — Gene {segment}")
    plt.xlabel(f"Frequência Relativa Mediana ({segment} Sem Outlier)")
    plt.ylabel(f"Frequência Relativa Mediana ({segment} Com Outlier)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    plt.show()


# ---------------------------------------------------------------------------- #
#                              EXPORTAÇÃO DE ARQUIVOS                          #
# ---------------------------------------------------------------------------- #

def salvar_fasta_filtrado(df: pd.DataFrame, output_path: str):
    with open(output_path, "w") as f:
        for _, row in df.iterrows():
            header = (f">{row['ID']} | Subtype:{row['Subtype']} | Species:{row['Species']} "
                      f"| Location:{row['Location']} | Data:{row['Data']} "
                      f"| Segment:{row['Segment']} | Class:{row['Class']}")
            f.write(header + "\n")
            f.write(row['Sequence'] + "\n")
    print(f"FASTA salvo em: {output_path}")


def salvar_fastas_combinados(df: pd.DataFrame, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    print("\n" + "=" * 50)
    print("Gerando arquivos FASTA combinados")
    print("=" * 50)
    for (subtype, segment), group in df.groupby(['Subtype', 'Segment']):
        if subtype == 'Unknown':
            continue
        file_name = f"{subtype}({segment})_todos.fa"
        file_path = os.path.join(output_dir, file_name)
        with open(file_path, "w") as f:
            for _, row in group.iterrows():
                f.write(f">{row['ID']}\n")
                seq = row['Sequence']
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        print(f"Criado: {file_name} ({len(group)} sequências)")
    print(f"Arquivos salvos em: {output_dir}")
    print("=" * 50)


# ---------------------------------------------------------------------------- #
#                                 PRINTS AUXILIARES                            #
# ---------------------------------------------------------------------------- #

def printar_secao(titulo: str):
    borda = "-" * 78
    print(f"\n{borda}")
    print(f"{titulo:^78}")
    print(borda)


def relatorio_remocoes(label: str, antes: int, depois: int):
    print(f"{label}: {antes - depois} removidas. Restantes: {depois}")



def main():
    data = carregar_dados(DOWNLOAD_PATH)
    df   = construir_dataframe(data)
    print(f"Total de sequências carregadas: {len(df)}")

    
    plotar_countplot(df, x='Class', titulo='Contagem de Sequências por Hospedeiro',
                     xlabel='Hospedeiro', ylabel='Número de Sequências')
    print("Contagem por hospedeiro:")
    print(df['Class'].value_counts().to_string())

    df["Data"] = df["Data"].apply(lambda x: x[-4:] if len(x) > 4 else x)
    plotar_countplot_por_segmento(df, ['HA', 'NA'])

    
    df.to_csv(os.path.join(DOWNLOAD_PATH, "csvcerto2.csv"), index=False)

    
    antes = len(df)
    df = corrigir_colunas(df)
    df = limpar_datas(df)
    relatorio_remocoes("Informações não padronizadas removidas", antes, len(df))

    plotar_countplot_por_segmento(df, ['HA', 'NA'])

    antes = len(df)
    df = remover_duplicatas(df)
    relatorio_remocoes("Duplicatas removidas", antes, len(df))

    antes = len(df)
    df = filtrar_location_invalida(df)
    relatorio_remocoes("Locations inválidas removidas", antes, len(df))

    
    pkl_path = os.path.join(HERE, "df_analise_N.pkl")
    try:
        df.to_pickle(pkl_path)
        print(f"PKL salvo em: {pkl_path}")
    except Exception as e:
        print(f"Erro ao salvar PKL: {e}")

    antes = len(df)
    df = filtrar_stop_codons(df)
    relatorio_remocoes("Sequências com stop codons inválidos removidas", antes, len(df))

    salvar_fasta_filtrado(df, os.path.join(DOWNLOAD_PATH, "df_sequencias_corretas.fa"))

    df = limpar_location(df)
    antes = len(df)
    df = filtrar_sequencias_invalidas(df)
    relatorio_remocoes("Sequências com caracteres inválidos removidas", antes, len(df))

    contagem_cidades = df["Location"].value_counts()
    print(contagem_cidades)
    print(f"Quantidade de locais únicos: {len(contagem_cidades)}")

    df = geocodificar_cidades(df, CACHE_FILE)
    antes = len(df)
    df.dropna(subset=['Country'], inplace=True)
    relatorio_remocoes("Linhas com país nulo removidas", antes, len(df))

    df = adicionar_continente(df)
    df['Data'] = pd.to_numeric(df['Data'], errors='coerce').astype(int)

    plotar_countplot(df, x='Subtype', titulo='Contagem de Sequências por Subtipo',
                     xlabel='Subtipo', ylabel='Número de Sequências', figsize=(15, 8))

    for hospedeiro in df['Class'].unique():
        subset = df[df['Class'] == hospedeiro]
        plotar_countplot(subset, x='Subtype',
                         titulo=f'Contagem de Subtipos — Hospedeiro: {hospedeiro}',
                         xlabel='Subtipo', ylabel='Número de Sequências', figsize=(10, 6))

    for seg in ["HA", "NA"]:
        plotar_cumulativo(df, seg, "Country", "Country",
                          f"Curva Cumulativa ({seg}) por Data e País")
        plotar_cumulativo(df, seg, "Continente", "Continente",
                          f"Curva Cumulativa ({seg}) por Data e Continente")

    plotar_cumulativo_geral(df)

    for seg in ["HA", "NA"]:
        plotar_cumulativo(df, seg, "Class", "Hospedeiro",
                          f"Curva Cumulativa ({seg}) por Data e Hospedeiro")

    df.drop(columns=["Continente"], inplace=True)
    print("Coluna 'Continente' removida do DataFrame final.")

    salvar_fastas_combinados(df, os.path.join(HERE, 'FastaCombinados'))

    df["Nucleotides_count"] = df["Sequence"].apply(len)
    df["Codons_count"]      = df["Sequence"].apply(lambda seq: len(seq) // 3)

    if not COM_CODONS:
        return

    df1 = df[df["Sequence"].str.fullmatch("[AGCT]+")].copy()

    for seg in ["HA", "NA"]:
        df_seg = df1[df1["Segment"] == seg]
        print(f'\nContagem por hospedeiro — gene {seg} (após todas as filtragens):')
        print(df_seg['Class'].value_counts().to_string())

    print(f"\nTotal de sequências para perfil de uso de códons: {len(df1)}")

    codon_df_raw, codon_df_norm = gerar_tabela_codons(df1, normalizar=True)

    plotar_countplot(df1, x='Subtype',
                     titulo='Contagem de Sequências por Subtipo (análise de códons)',
                     xlabel='Subtipo', ylabel='Número de Sequências', figsize=(15, 8))

    relatorio_stop_codons(df1, codon_df_norm)

    df_codons = pd.concat([df1.reset_index(drop=True), codon_df_norm], axis=1)

    for seg in ["NA", "HA"]:
        codon_dfs_seg = df_codons[df_codons["Segment"] == seg]
        plotar_boxplot_codons(codon_dfs_seg, codon_df_norm, seg, DOWNLOAD_PATH)

    codon_usage_table = pd.concat([df1["ID"].reset_index(drop=True), codon_df_norm], axis=1)
    codon_usage_table.to_csv(os.path.join(DOWNLOAD_PATH, "codon_usage_table.csv"), index=False)
    print(f"Tabela de uso de códons salva.")

    Tabela_De_Outliers, df_com_out, df_sem_out, outlier_mask = detectar_outliers_por_segmento(
    df1.reset_index(drop=True), codon_df_norm)

    Tabela_De_Outliers.to_csv(os.path.join(DOWNLOAD_PATH, "Tabela_De_Outliers.csv"), index=False)
    df_com_out.to_csv(os.path.join(DOWNLOAD_PATH, "df_com_out.csv"), index=False)
    df_sem_out.to_csv(os.path.join(DOWNLOAD_PATH, "df_sem_out.csv"), index=False)

    printar_secao("DataFrame básico")
    print(df)
    printar_secao("df1 — apenas AGCT")
    print(df1)
    printar_secao("Tabela de Outliers")
    print(Tabela_De_Outliers)
    printar_secao("Sequências com outliers")
    print(df_com_out)
    df_com_out.info()
    printar_secao("Sequências sem outliers")
    print(df_sem_out)
    df_sem_out.info()

    todos_com = {}
    todos_sem = {}

    for seg in ["HA", "NA"]:
        com = codons_por_grupo(df1.reset_index(drop=True), codon_df_raw, df_com_out["ID"], seg)
        sem = codons_por_grupo(df1.reset_index(drop=True), codon_df_raw, df_sem_out["ID"], seg)
        todos_com[seg] = com
        todos_sem[seg] = sem

        printar_secao(f"Todos os códons — COM outliers ({seg})")
        print(com)
        com.info()

        printar_secao(f"Boxplot — COM outliers ({seg})")
        plotar_boxplot_simples(
            com.drop(columns=["ID", "Segment"], errors='ignore'),
            titulo=f"Box Plot — Influenza A ({seg}) com Outliers"
        )

        printar_secao(f"Todos os códons — SEM outliers ({seg})")
        print(sem)

        printar_secao(f"Boxplot — SEM outliers ({seg})")
        plotar_boxplot_simples(
            sem.drop(columns=["ID", "Segment"], errors='ignore'),
            titulo=f"Box Plot — Influenza A ({seg}) sem Outliers"
        )

        com_vals = com.drop(columns=["ID", "Segment"], errors='ignore')
        sem_vals = sem.drop(columns=["ID", "Segment"], errors='ignore')
        comparacao = pd.DataFrame({"Não outlier": sem_vals.mean(), "Outlier": com_vals.mean()})

        printar_secao(f"Comparação com e sem outliers — {seg}")
        plotar_comparacao_barras(comparacao, seg)
        plotar_dispersao_medias(comparacao, seg)

    todos_com_df = pd.concat(todos_com.values())
    todos_sem_df = pd.concat(todos_sem.values())
    plotar_correlacao_ha_na(todos_sem_df, todos_com_df)

    for seg in ["HA", "NA"]:
        com_vals = todos_com[seg].drop(columns=["ID", "Segment"], errors='ignore')
        sem_vals = todos_sem[seg].drop(columns=["ID", "Segment"], errors='ignore')
        plotar_impacto_outliers(sem_vals.median(), com_vals.median(), seg)

    printar_secao("Gerando FASTA — sequências com outliers")
    seq_com_out = df1[df1["ID"].isin(df_com_out["ID"])]
    fasta_out   = os.path.join(DOWNLOAD_PATH, "filtrado.fa")
    with open(fasta_out, "w") as f:
        for _, row in seq_com_out.iterrows():
            f.write(f">{row['ID']}\n{row['Sequence']}\n")
    print(f"FASTA salvo em: {fasta_out}")


if __name__ == "__main__":
    main()