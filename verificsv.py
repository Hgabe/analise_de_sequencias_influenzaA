import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import matplotlib.patches as mpatches
from geopy.geocoders import Nominatim
import time
import pycountry_convert as pc
import json
import pycountry
from geopy.exc import GeocoderTimedOut, GeocoderServiceError
from googletrans import Translator
here = os.path.abspath(os.path.dirname(__file__))

download_path = os.path.join(here, 'ArqFasta')

comcodons = True
data = []


identifier_pattern = re.compile(
    r"(?P<ID>[^\s]+)\s(?P<Type>[A-Z])\/(?P<Species>[^/]+(?:\s[^/]+)*)\/(?P<Location>[^\/]+)\/(?P<Data>[^/]+)\/(?P<Sequence_type>[^\s]+)\s"
)
()
for file_name in os.listdir(download_path):
    if file_name.endswith(".fa"):
        file_path = os.path.join(download_path, file_name)

        # Captura do nome da classe assim: (x)
        class_match = re.search(r"\((.*?)\)", file_name)
        class_name = class_match.group(1) if class_match else "Unknown"

        # Captura do nome do segmento assim: _x_
        segment_match = re.search(r"_(.*?)_", file_name)
        segment = segment_match.group(1) if segment_match else "Unknown"

        # Captura o subtipo no nome
        Subtype_match = re.search(r"^(H\d+N\d+)", file_name)
        Subtype = Subtype_match.group(1) if Subtype_match else "Unknown"

        
        with open(file_path, "r") as fasta_file:
            identifier = ""
            sequence = ""

            for line in fasta_file:
                line = line.strip()
                if line.startswith(">"):
                    if identifier and sequence:
                        data.append([identifier, sequence, class_name, segment, Subtype])
                    identifier = line[1:]
                    sequence = ""
                else:
                    sequence += line

            if identifier and sequence:
                data.append([identifier, sequence, class_name, segment, Subtype])


data_cleaned = []
for entry in data:
    match = identifier_pattern.match(entry[0])
    if match:
        ID = match.group("ID")
        Type = match.group("Type")
        Species = match.group("Species")
        Location = match.group("Location")
        Data = match.group("Data")
        Sequence_type = match.group("Sequence_type")

        data_cleaned.append([
            ID, Type, Species, Location, Data, Sequence_type,
            entry[1], entry[2], entry[3], entry[4]
        ])


df = pd.DataFrame(data_cleaned, columns=[
    "ID", "Type", "Species", "Location", "Sequence_type", "Data",
    "Sequence", "Class", "Segment", "Subtype"
])
print(f"Total de sequências baixadas: {len(df)}")

print(f'Contagem de Sequências por hospedeiro')
for hospedeiro, count in df['Class'].value_counts().items():
    print(f'{hospedeiro}: {count}')
plt.figure(figsize=(6, 4))
sns.countplot(
    x='Class',
    data=df,
    color='gray',
    order=df['Class'].value_counts().index
)
plt.title(f'Contagem de Sequências por Hospedeiro', fontsize=16)
plt.xlabel('Hospedeiro', fontsize=12)
plt.ylabel('Número de Sequências', fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.yscale("log")
plt.xticks(rotation=90)
plt.tight_layout() 
plt.show()
atual = len(df)
df["Data"] = df["Data"].apply(lambda x: x[-4:] if len(x) > 4 else x)
segmentos_alvo = ['HA', 'NA']

for segmento in segmentos_alvo:
    df_segmento = df[df['Segment'] == segmento]
    contagens = df_segmento['Class'].value_counts()
    print(f'\nContagem de Sequências por Hospedeiro para o gene {segmento}')
    print(contagens)
    plt.figure(figsize=(8, 5))
    sns.countplot(
        x='Class',
        color='gray',
        data=df_segmento,
        order=contagens.index # Ordena as barras da mais alta para a mais baixa
    )
    
    plt.title(f'Contagem de Sequências por Hospedeiro para o gene {segmento}', fontsize=14)
    plt.xlabel('Hospedeiro', fontsize=12)
    plt.ylabel('Número de Sequências', fontsize=12)
    plt.yscale("log") 
    plt.xticks(rotation=90)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()
csv_path = os.path.join(download_path, "csvcerto2.csv")
df.to_csv(csv_path, index=False)
print(f"Arquivo CSV salvo em: {csv_path}")


df[["Data", "Sequence_type", "Location", "Species"]] = df.apply(
    lambda row: pd.Series([
        row["Sequence_type"] if "/" in row["Data"] else row["Data"],
        row["Location"] if "/" in row["Data"] else row["Sequence_type"],
        row["Species"] if "/" in row["Data"] else row["Location"],
        np.nan if "/" in row["Data"] else row["Species"]
    ]),
    axis=1
)


df["Data"] = df["Data"].apply(lambda x: x[-4:] if len(x) > 4 else x)
df["Data"] = df["Data"].str.replace(r"^([0-9]{2})$", r"19\1", regex=True)
df = df[df["Data"].str.len() != 1]
df = df[df["Data"].str.len() != 3]
df["Data"] = df["Data"].apply(lambda x: x[1:] if len(x) == 5 else x)
df["Data"] = df["Data"].apply(lambda x: x if x.isdigit() else None)

df = df.dropna(subset=["Data"])
#aqui todas as linhas swine NA foram removidas
print(f"Linhas com informações não padronizadas removidas. Total de sequências removidas: {atual-len(df)}")
print (f"Total de sequências restantes: {len(df)}")
atual = len(df)
df = df.drop_duplicates(subset=["ID", "Sequence"], keep='last')
segmentos_alvo = ['HA', 'NA']

for segmento in segmentos_alvo:
    df_segmento = df[df['Segment'] == segmento]
    contagens = df_segmento['Class'].value_counts()
    print(f'\nContagem de Sequências por Hospedeiro para o gene {segmento}')
    print(contagens)
    plt.figure(figsize=(8, 5))
    sns.countplot(
        x='Class',
        color='gray',
        data=df_segmento,
        order=contagens.index # Ordena as barras da mais alta para a mais baixa
    )
    
    plt.title(f'Contagem de Sequências por Hospedeiro para o gene {segmento}', fontsize=14)
    plt.xlabel('Hospedeiro', fontsize=12)
    plt.ylabel('Número de Sequências', fontsize=12)
    plt.yscale("log") 
    plt.xticks(rotation=90)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()
print(f"Linhas duplicadas removidas. Total de sequências removidas: {atual-len(df)}")
print (f"Total de sequências restantes: {len(df)}")

atual = len(df)
df = df[~(df["Location"].str.contains(r"\d", regex=True) | (df["Location"].str.len() <= 2))]
print(f"linhas removidas por location inválida: {atual-len(df)}")
print(f"Total de sequências restantes: {len(df)}")

pkl_output_name = "df_analise_N.pkl"
pkl_output_path = os.path.join(here, pkl_output_name)

try:
    df.to_pickle(pkl_output_path)
    print(f"\nDataFrame final (df) salvo como arquivo PKL em: {pkl_output_path}")
except Exception as e:
    print(f"Erro ao salvar o arquivo PKL: {e}")
def count_stop_codons(seq):
    count = 0
    stop_codons = ["TAA", "TGA", "TAG"]
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in stop_codons:
            count += 1
    return count

# Crie a coluna com a contagem de códons de parada
df["StopCodons_count"] = df["Sequence"].apply(count_stop_codons)

# Defina a condição principal que todas as sequências devem atender
mask_comeca_com_atg = df['Sequence'].str.startswith('ATG')

# Defina as duas lógicas para a parte dos códons de parada
mask_sem_stop_codon = (df['StopCodons_count'] == 0)
mask_com_stop_codon = (df['StopCodons_count'] == 1) & \
                      (df['Sequence'].str.endswith(('TAA', 'TAG', 'TGA')))

# Combine todas as máscaras:
# A sequência precisa COMEÇAR com 'ATG' E (ter 0 códons de parada OU ter 1 códon e terminar com um deles)
mask_final = mask_comeca_com_atg & (mask_sem_stop_codon | mask_com_stop_codon)
atual = len(df)
# Aplique a máscara final ao DataFrame
df = df[mask_final].copy()

# Remova a coluna temporária "StopCodons_count"
df.drop(columns=["StopCodons_count"], inplace=True)

print(f"Linhas com stop codons em posições incorretas, iniciadas incorretamente ou com mais de um stop codon removidas: {atual-len(df)}")
print (f"Total de sequências restantes: {len(df)}")


# --- Código para salvar o DataFrame como arquivo FASTA ---

# Define o nome do arquivo FASTA de saída
fasta_output_name = "df_sequencias_corretas.fa"
fasta_output_path = os.path.join(download_path, fasta_output_name)

# Abre o arquivo para escrita
with open(fasta_output_path, "w") as fasta_file:
    # Itera sobre cada linha (sequência) no DataFrame processado
    for index, row in df.iterrows():
        # Constrói o cabeçalho FASTA (linha que começa com >)
        # Incluindo informações relevantes para identificação.
        header = f">{row['ID']} | Subtype:{row['Subtype']} | Species:{row['Species']} | Location:{row['Location']} | Data:{row['Data']} | Segment:{row['Segment']} | Class:{row['Class']}"
        
        # Escreve o cabeçalho no arquivo
        fasta_file.write(header + "\n")
        
        # Escreve a sequência no arquivo
        # É comum que sequências sejam quebradas em blocos de 60 ou 80 caracteres em arquivos FASTA,
        # mas para simplicidade, vamos escrever a sequência inteira em uma linha.
        fasta_file.write(row['Sequence'] + "\n")

print(f"\nDataFrame salvo como arquivo FASTA em: {fasta_output_path}")

# --- Fim do código ---

df['Location'] = df['Location'].str.replace(r'-\s?[Cc]$', '', regex=True, flags=re.IGNORECASE)
df['Location'] = df['Location'].str.replace(r'-\s?SB$', '', regex=True, flags=re.IGNORECASE)
df['Location'] = df['Location'].str.strip()
total_antes = len(df)
def percentual_invalidos(seq):
    total = len(seq)
    invalidos = sum(1 for c in seq if c not in "AGCT")
    return invalidos / total if total > 0 else 1 

df["Invalid_pct"] = df["Sequence"].apply(percentual_invalidos)

df = df[df["Invalid_pct"] < 0.01]
total_depois = len(df)
total_real = total_antes - total_depois
print(f"Número total de sequências removidas por caracteres inválidos: {total_real}")

contagem_cidades = df["Location"].value_counts()
print(contagem_cidades)
total_De_cidades = len(contagem_cidades)
print(f"Quantidade de nomes de cidades que aparecem: {total_De_cidades}")



cache_file = "cache_cidade_pais.json"
geolocator = Nominatim(user_agent="city-to-country-script-v5-final")
translator = Translator()
traducao_cache = {} 

cidade_para_pais_dict = {}
if os.path.exists(cache_file) and os.path.getsize(cache_file) > 0:
    try:
        with open(cache_file, "r", encoding="utf-8") as f:
            cidade_para_pais_dict = json.load(f)
    except Exception as e:
        print("Erro ao carregar cache, recriando:", e)
        cidade_para_pais_dict = {}

def get_country_name(country_code):
    try:
        country = pycountry.countries.get(alpha_2=country_code)
        if country:
            return country.name
    except Exception:
        return None
    return None

def traduzir_nome(nome):
    if nome in traducao_cache:
        return traducao_cache[nome]
    try:
        translated = translator.translate(nome, dest='en').text
        traducao_cache[nome] = translated
        print(f"Traduzido '{nome}' -> '{translated}'")
        return translated
    except Exception:
        return nome


def geocode_with_nominatim(cidade, use_translation=False, timeout=10):
    if use_translation:
        cidade = traduzir_nome(cidade)
        
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


mock_llm_responses = {
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
    "Osaka": "Japan","Osaka-c": "Japan", "Pingtung": "Taiwan", "Raesefeld": "Germany", "Shaanxi": "China",
    "Sichuan": "China", "SouthDakota": "United States", 'an': None, "Taichung": "Taiwan", "Tianjin": "China",
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
    "Tianjinnankai": "China", "Tianjinninghe": "China", "Tianjintanggu": "China", "Tientsin": "China",
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
    "Shijiazhuang": "China", "Taipei City": "Taiwan", "Yunlin": "Taiwan", "Jiang Xi": "China",
    "Nantou": "Taiwan", "Yilan": "Taiwan", "Djankoy": "Ukraine", "bird feces": None, "Taishun": "China",
    "Khunt lake": "Mongolia", "Sanmenxia": "China", "kumamoto": "Japan", "Egypt ": "Egypt", "Yuhang": "China",
    "Nanji": "China", "San Jiang": "China", "Hei Longjiang": "China", "Tumuji": "China", "Antarctic": "Antarctica",
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
    "GuangzhouSB": "China", "Azarbayejan-Sharghi":"Iran"
}


def cidade_para_pais_master(cidade, timeout=10):
    """
    Função principal com lógica de fallback aprimorada.
    """

    pais_nominatim = geocode_with_nominatim(cidade, use_translation=False, timeout=timeout)
    if pais_nominatim:
        print(f"Geocodificado diretamente com Nominatim: '{cidade}' -> '{pais_nominatim}'")
        return pais_nominatim
  
    pais_nominatim_traduzido = geocode_with_nominatim(cidade, use_translation=True, timeout=timeout)
    if pais_nominatim_traduzido:
        print(f"Geocodificado com tradução: '{cidade}' -> '{pais_nominatim_traduzido}'")
        return pais_nominatim_traduzido
        
   
    if cidade in mock_llm_responses:
        pais_llm = mock_llm_responses[cidade]
        if pais_llm:
            print(f"LLM corrigiu '{cidade}' -> '{pais_llm}'")
            return pais_llm
    
    print(f"Nenhum resultado para '{cidade}'")
    return None



cidades_unicas = df["Location"].dropna().unique()
to_check = [c for c in cidades_unicas if c not in cidade_para_pais_dict or cidade_para_pais_dict[c] is None]

max_tries = 3
save_every = 50

for i, cidade in enumerate(to_check, 1):
    if cidade not in cidade_para_pais_dict:
        cidade_para_pais_dict[cidade] = None

    if cidade_para_pais_dict[cidade] is not None:
        continue

    cidade_str = str(cidade)
    if cidade_str.isdigit() or len(cidade_str) <= 2:
        cidade_para_pais_dict[cidade] = None
        continue

    print(f"[{i}/{len(to_check)}] Consultando: {cidade_str}")

    success = False
    for attempt in range(1, max_tries + 1):
        try:
            pais = cidade_para_pais_master(cidade_str)
            cidade_para_pais_dict[cidade] = pais
            success = True
            break
        except GeocoderTimedOut:
            wait = 2 ** attempt
            print(f"Timeout na tentativa {attempt}. Backoff {wait}s...")
            time.sleep(wait)
        except GeocoderServiceError as e:
            print(f"Erro de serviço: {e}. Esperando 10s antes de tentar novamente...")
            time.sleep(10)
        except Exception as e:
            print(f"Erro inesperado consultando '{cidade_str}': {e}")
            break

    if not success and cidade_para_pais_dict.get(cidade) is None:
        cidade_para_pais_dict[cidade] = None

    time.sleep(1)

    if i % save_every == 0:
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(cidade_para_pais_dict, f, ensure_ascii=False, indent=2)
        print(f"Salvo progresso ({i} consultas).")


with open(cache_file, "w", encoding="utf-8") as f:
    json.dump(cidade_para_pais_dict, f, ensure_ascii=False, indent=2)
print("Cache atualizado e salvo.")

df["Country"] = df["Location"].map(cidade_para_pais_dict)
atual = len(df)
df.dropna(subset=['Country'], inplace=True)
print(f"Linhas com país nulo removidas. Total de sequências removidas: {atual-len(df)}")
print (f"Total de sequências restantes: {len(df)}")

continente_cache = {}

def pais_para_continente(pais):
    if pais in continente_cache:
        return continente_cache[pais]
    try:
        country_code = pc.country_name_to_country_alpha2(pais, cn_name_format="default")
        continent_code = pc.country_alpha2_to_continent_code(country_code)
        continent_name = {
            'NA': 'América do Norte',
            'SA': 'América do Sul',
            'AS': 'Ásia',
            'EU': 'Europa',
            'AF': 'África',
            'OC': 'Oceania',
            'AN': 'Antártica'
        }.get(continent_code, 'Outros')
        
        continente_cache[pais] = continent_name
        return continent_name
    except Exception:
        continente_cache[pais] = 'Outros'
        return 'Outros'


df["Continente"] = df["Country"].apply(pais_para_continente)
df['Data'] = pd.to_numeric(df['Data'])
df['Data'] = df['Data'].astype(int)


plt.figure(figsize=(15, 8))
sns.countplot(
    x='Subtype',
    data=df,
    color='gray',
    order=df['Subtype'].value_counts().index
)
plt.title('Contagem de Sequências por Subtipo', fontsize=16)
plt.xlabel('Subtipo', fontsize=12)
plt.ylabel('Número de Sequências', fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.yscale("log")
plt.xticks(rotation=90)
plt.tight_layout() 
plt.show()


for Hospedeiro in df['Class'].unique():
    plt.figure(figsize=(10, 6))
    subset = df[df['Class'] == Hospedeiro]
    sns.countplot(
        x='Subtype',
        data=subset,
        color='gray',
        order=subset['Subtype'].value_counts().index
    )
    plt.title(f'Contagem de subtipos por Subtipo para Hospedeiro(teste): {Hospedeiro}', fontsize=16)
    plt.xlabel('Subtipo', fontsize=12)
    plt.ylabel('Número de Sequências', fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.yscale("log")
    plt.xticks(rotation=90)
    plt.tight_layout() 
    plt.show()

# ("------------------------------------------------------------------------------")
# ("                            gerando gráficos cumulativos                      ")
# ("                               de sequencias NA e HA                          ")
# ("                                  por data e país                             ")
# ("------------------------------------------------------------------------------")
for segmento in ["HA", "NA"]:
    df_filtrado = df[df["Segment"] == segmento]
    df_counts = df_filtrado.groupby(["Data", "Country"]).size().unstack(fill_value=0)
    
    df_counts.index.name = "Data de Coleta"
    df_counts.columns.name = "Country"
    df_cumsum = df_counts.cumsum()
    df_cumsum = df_cumsum.apply(pd.to_numeric, errors="coerce").fillna(0)
    df_cumsum.plot(
        kind="line",
        stacked=False,
        title=f"Curva Cumulativa de Sequências ({segmento}) por Data de Coleta e País"
    )
    plt.xticks(
        ticks=df_cumsum.index,
        labels=df_cumsum.index,
        rotation=90,
        ha='center',
        fontsize=7
    )
    plt.yscale("log")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel("Data de Coleta")
    plt.ylabel("Quantidade Acumulada (log)")
    plt.tight_layout()
    plt.show()

# ("------------------------------------------------------------------------------")
# ("                            Gerando gráficos cumulativos                      ")
# ("                               de sequencias NA e HA e                        ")
# ("                                   por continente                             ")
# ("------------------------------------------------------------------------------")
for segmento in ["HA", "NA"]:
    df_filtrado_cont = df[df["Segment"] == segmento]
    df_counts_cont = df_filtrado_cont.groupby(["Data", "Continente"]).size().unstack(fill_value=0)
    
    df_counts_cont.index.name = "Data de Coleta"
    df_counts_cont.columns.name = "Continente"
    df_cumsum_cont = df_counts_cont.cumsum()
    df_cumsum_cont = df_cumsum_cont.apply(pd.to_numeric, errors="coerce").fillna(0)
    df_cumsum_cont.plot(
        kind="line",
        stacked=False,
        title=f"Curva Cumulativa de Sequências do vírus Influenza A para o gene {segmento} por Continente"
    )
    all_dates = df_cumsum.index.unique()
    plt.xticks(
        ticks=df_cumsum.index,
        labels=df_cumsum.index,
        rotation=90,
        ha='center',
        fontsize=7
    )
    plt.yscale("log")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel("Data de Coleta")
    plt.ylabel("Quantidade Acumulada (log)")
    plt.tight_layout()
    plt.show()

# A coluna Continente foi criada para análises intermediárias, mas foi removida do DataFrame final
df.drop(columns=["Continente"], inplace=True)
print("Coluna 'Continente' removida do DataFrame final.")


# --- Bloco de Geração de FASTA Combinado ---
# O caminho para salvar os novos arquivos FASTA
output_fasta_path = os.path.join(here, 'FastaCombinados')
os.makedirs(output_fasta_path, exist_ok=True) # Cria o diretório se não existir

print("\n" + "="*50)
print("Iniciando a geração de arquivos FASTA combinados")
print("="*50)

# Agrupando o DataFrame por 'Subtype' e 'Segment'
# Ex: (H1N1, HA), (H1N1, NA), (H5N1, HA), etc.
grouped = df.groupby(['Subtype', 'Segment'])

for (subtype, segment), group in grouped:
    # Ignora 'Subtype's que são 'Unknown' (você pode ajustar isso)
    if subtype == 'Unknown':
        continue
        
    # Nome do novo arquivo FASTA
    # Ex: H1N1(HA)_todos.fa ou H5N1(NA)_todos.fa
    output_file_name = f"{subtype}({segment})_todos.fa"
    output_file_path = os.path.join(output_fasta_path, output_file_name)
    
    # Abrindo o arquivo para escrita
    with open(output_file_path, "w") as outfile:
        # Iterando sobre as linhas do grupo (que contém apenas as sequências
        # do Subtype e Segmento atual)
        for index, row in group.iterrows():
            # O cabeçalho FASTA é o 'ID' (que é a linha original que começa com '>').
            # O ID completo é mais informativo do que apenas usar o Subtype.
            fasta_header = f">{row['ID']}"
            
            # A sequência
            sequence = row['Sequence']
            
            # Escreve o cabeçalho e a sequência no novo arquivo
            outfile.write(fasta_header + "\n")
            
            # Escreve a sequência, garantindo quebras de linha a cada 60-80 caracteres
            # para o formato FASTA padrão (opcional, mas boa prática)
            for i in range(0, len(sequence), 80):
                outfile.write(sequence[i:i+80] + "\n")
                
    print(f"Criado arquivo: {output_file_name} com {len(group)} sequências.")

print("="*50)
print("Geração de arquivos FASTA concluída.")
print(f"Arquivos salvos em: {output_fasta_path}")
print("="*50)






# ("------------------------------------------------------------------------------")
# ("                            Gerando gráfico cumulativo                        ")
# ("                             por sequencias NA e HA e                         ")
# ("                                     por data                                 ")
# ("------------------------------------------------------------------------------")
df_counts = df.groupby(["Data", "Segment"]).size().unstack(fill_value=0)
df_counts.index.name = "Data de Coleta"
df_counts.columns.name = "Segmento"
df_cumsum = df_counts.cumsum()
df_cumsum.plot(kind="line", stacked=False, color=["blue", "red"], title="Curva Cumulativa de Sequências por Data de Coleta e Segmento")
plt.xticks(
    ticks=df_cumsum.index,
    labels=df_cumsum.index,
    rotation=90,
    ha='center',
    fontsize=7
)
plt.yscale("log")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel("Data de Coleta")
plt.ylabel("Quantidade Acumulada (log)")
plt.tight_layout()
plt.show()

# ("------------------------------------------------------------------------------")
# ("                          Gerando gráficos cumulativos                        ")
# ("                             de sequencias NA e HA                            ")
# ("                           por hospedeiro e por data                          ")
# ("------------------------------------------------------------------------------")

for segmento in ["HA", "NA"]:
    df_filtrado_seg = df[df["Segment"] == segmento]
    df_counts = df_filtrado_seg.groupby(["Data", "Class"]).size().unstack(fill_value=0)
    df_counts.index.name = "Data de Coleta"
    df_counts.columns.name = "Hospedeiro"
    df_cumsum = df_counts.cumsum()
    df_cumsum.plot(kind="line", stacked=False, color=["blue", "red", "green"], title=f"Curva Cumulativa de Sequências do vírus Influenza A para o gene {segmento} por Data de Coleta e Hospedeiro")
    plt.xticks(
        ticks=df_cumsum.index,
        labels=df_cumsum.index,
        rotation=90,
        ha='center',
        fontsize=7
    )
    plt.yscale("log")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel("Data de Coleta")
    plt.ylabel("Quantidade Acumulada (log)")
    plt.tight_layout()
    plt.show()




df["Nucleotides_count"] = df["Sequence"].apply(lambda seq: len (seq))
df["Codons_count"] = df["Sequence"].apply(lambda seq: len (seq)//3)
if comcodons:
    def contar_codons(seq):
        codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
        return dict(Counter(codons))



    df1 = df[df["Sequence"].str.fullmatch("[AGCT]+")]
    for segmento in ["HA", "NA"]:
        df_segmento = df1[df1["Segment"] == segmento]
        print(f'Contagem de Sequências por hospedeiro para o gene {segmento} após todas as filtragens:')
        for hospedeiro, count in df_segmento['Class'].value_counts().items():
            print(f'{hospedeiro}: {count}')
            
        codon_dicts = df1["Sequence"].apply(contar_codons)
        print(f"total de sequencias usadas para perfil de uso de codons: {len(df1)}" )
        stop_codons = ["TAA", "TGA", "TAG"]

    def count_stop_codons(seq):
        count = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in stop_codons:
                count += 1
        return count
    
    plt.figure(figsize=(15, 8))
    sns.countplot(
        x='Subtype',
        data=df1,
        color='gray',
        order=df1['Subtype'].value_counts().index
    )
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.yscale("log")
    plt.title('Contagem de Sequências por Subtipo para análise de uso de códon', fontsize=16)
    plt.xlabel('Subtipo', fontsize=12)
    plt.ylabel('Número de Sequências', fontsize=12)

   
    plt.xticks(rotation=90, ha='right')

    
    plt.tight_layout()  
    plt.show()

    df1_copy = df1.copy()
    df1_copy["StopCodons_count"] = df1_copy["Sequence"].apply(count_stop_codons)
    zero_stop = len(df1_copy[df1_copy["StopCodons_count"] == 0])
    one_stop_codon_count = len(df1_copy[df1_copy["StopCodons_count"] == 1])
    more_than_one_stop_codon_counts = df1_copy[df1_copy["StopCodons_count"] > 1]["StopCodons_count"].value_counts().sort_index()
    quantidade_total_com_mais = sum(more_than_one_stop_codon_counts)
    print(f"Número de sequências com 0 códons de parada: {zero_stop}")
    print(f"Número de sequências com exatamente 1 codon de parada: {one_stop_codon_count}")
    print(f"\nContagem de sequências com mais de 1 codon de parada:{quantidade_total_com_mais}")
    print("\nContagem de sequências com mais de 1 codon de parada para cada quantidade:")
    print(more_than_one_stop_codon_counts.to_string())

    codon_df = pd.DataFrame.from_records(codon_dicts).fillna(0).astype(int)
    codon_df.drop(columns=["TAA", "TGA", "TAG"], errors='ignore', inplace=True)
    mask_one_stop_codon = (df1_copy["StopCodons_count"] == 1)

    # Filtre o DataFrame usando a máscara
    df_one_stop = df1_copy[mask_one_stop_codon].copy()

    # Crie uma máscara para cada códon de parada no final da sequência.
    mask_taa_end = df_one_stop["Sequence"].str.endswith('TAA')
    mask_tag_end = df_one_stop["Sequence"].str.endswith('TAG')
    mask_tga_end = df_one_stop["Sequence"].str.endswith('TGA')

    # Conte o número de sequências que terminam com cada códon.
    taa_count = len(df_one_stop[mask_taa_end])
    tag_count = len(df_one_stop[mask_tag_end])
    tga_count = len(df_one_stop[mask_tga_end])

    # Imprima os resultados
    print(f"Sequências terminadas em TAA: {taa_count}")
    print(f"Sequências terminadas em TAG: {tag_count}")
    print(f"Sequências terminadas em TGA: {tga_count}")
    codons_counts_aligned = df1["Codons_count"].reindex(codon_df.index)
    normalized_codon_df = codon_df.div(codons_counts_aligned, axis=0)

    # O restante do código usa a nova variável normalizada
    df_codons = pd.concat([df1.reset_index(drop=True), normalized_codon_df.reset_index(drop=True)], axis=1)
    codon_dfs = pd.concat([df1["Segment"].reset_index(drop=True), normalized_codon_df.reset_index(drop=True)], axis=1)
    
    def geraboxplot():
        estatisticas_codons = pd.DataFrame({
            "total": codon_df.sum(),
            "Media": codon_df.mean(),
            "DesvioPadrao": codon_df.std()
        })
        estatisticas_codons.to_csv(os.path.join(download_path, f"estatisticas_codons_{segment}.csv"))

        codon_df_sorted = codon_dfs[sorted(codon_df.columns)]
        sns.boxplot(data=codon_df_sorted, color="lightgray", boxprops=dict(facecolor='gray'), showfliers=False) 
        
        mean_values = codon_df_sorted.mean()
        std_values = codon_df_sorted.std()

        for i, (mean, std) in enumerate(zip(mean_values, std_values)):
            plt.text(i, std, '__', ha='center', va='center', color='blue', fontsize=10)
            plt.text(i, mean, '__', ha='center', va='center', color='red', fontsize=10)

        plt.title("Box Plot com Média e Desvio Padrão dos Códons de Sequências do Vírus Influenza A Para o Gene " + segment)
        plt.xlabel("Códons")
        plt.ylabel("Frequência Relativa")
        plt.xticks(rotation=90)

        media_patch = mpatches.Patch(color='red', label='Média')
        dp_patch = mpatches.Patch(color='blue', label='Desvio Padrão')
        mediana_patch = mpatches.Patch(color='black', label='mediana')
        plt.legend(handles=[media_patch, dp_patch, mediana_patch], loc='upper right', bbox_to_anchor=(1, 1))

        plt.tight_layout()
        plt.show()
    
    segment = "NA"
    codon_dfs = df_codons[df_codons["Segment"] == segment]
    geraboxplot()

    segment = "HA"
    codon_dfs = df_codons[df_codons["Segment"] == segment]
    geraboxplot()

    codon_usage_table = pd.concat([df1["ID"].reset_index(drop=True), normalized_codon_df.reset_index(drop=True)], axis=1)
    codon_csv_path = os.path.join(download_path, "codon_usage_table.csv")
    codon_usage_table.to_csv(codon_csv_path, index=False)
    std_series = df_codons.drop(columns=df.columns).std()

    print(f"Arquivo com ID e contagens de códons salvo em: {codon_csv_path}")
    
print ("------------------------------------------------------------------------------")
print ("                               Exibindo df basico                             ")
print ("------------------------------------------------------------------------------")
print(df)
print ("------------------------------------------------------------------------------")
print ("                          Exibindo df1 sem caracteres especiais               ")
print ("------------------------------------------------------------------------------")
print(df1)

def detectar_outliers_por_grupo(grupo):
    Q1 = grupo.quantile(0.25)
    Q3 = grupo.quantile(0.75)
    IQR = Q3 - Q1
    outlier_mask = (grupo < (Q1 - 1.5 * IQR)) | (grupo > (Q3 + 1.5 * IQR))
    return outlier_mask

# A máscara de outlier agora é aplicada no DataFrame normalizado
outlier_mask = normalized_codon_df.groupby(df1["Segment"]).transform(detectar_outliers_por_grupo)

outlier_df = normalized_codon_df[outlier_mask]
Tabela_De_Outliers = pd.concat([df1["ID"].reset_index(drop=True), outlier_df.reset_index(drop=True)], axis=1)

def formatar_valores(row):
    return [f"{val}({col})" for col, val in row.items() if pd.notna(val)]

linhas_formatadas = Tabela_De_Outliers.apply(formatar_valores, axis=1)

Tabela_De_Outliers = pd.DataFrame({
    "Valores": [", ".join(linha) for linha in linhas_formatadas]
})
Tabela_De_Outliers = Tabela_De_Outliers[Tabela_De_Outliers["Valores"].str.count(r"\(") >= 2]

Tabela_De_Outliers.to_csv(os.path.join(download_path, "Tabela_De_Outliers.csv"), index=False)



df_out_temp = normalized_codon_df[outlier_mask]

print ("------------------------------------------------------------------------------")
print ("                           Exibindo tabela out                               ")
print ("------------------------------------------------------------------------------")

print(Tabela_De_Outliers)

print ("------------------------------------------------------------------------------")
print ("                       Exibindo sequencias com outliers                       ")
print ("------------------------------------------------------------------------------")

df_out_temp = pd.concat([df1["ID"].reset_index(drop=True), df_out_temp.reset_index(drop=True)], axis=1)
df_out_temp["tem_out"] = (~df_out_temp.drop(columns="ID").isna().all(axis=1)).astype(int)
df_com_out = df_out_temp[df_out_temp['tem_out'] == 1].drop(columns='tem_out')

df_com_out.to_csv(os.path.join(download_path, "df_com_out.csv"), index=False)

print(df_com_out)
df_com_out.info()

print ("------------------------------------------------------------------------------")
print ("                       Exibindo sequencias sem outliers                       ")
print ("------------------------------------------------------------------------------")

df_sem_out = df_out_temp[df_out_temp['tem_out'] == 0].drop(columns='tem_out')

df_sem_out.to_csv(os.path.join(download_path, "df_sem_out.csv"), index=False)

print(df_sem_out)
df_sem_out.info()

print ("------------------------------------------------------------------------------")
print ("                              Exibindo todos os codons                        ")
print ("                             de sequencias com outliers                       ")
print ("------------------------------------------------------------------------------")
df_filtrado = df1[df1["ID"].isin(df_com_out["ID"])].drop_duplicates(subset=["ID"]).reset_index(drop=True)
todos_codons_com_out = pd.concat([df_filtrado[["ID", "Segment"]], codon_df.reset_index(drop=True)], axis=1)
todos_codons_com_out = todos_codons_com_out.dropna(subset=["ID"])
print(todos_codons_com_out)
todos_codons_com_out.info()

print ("------------------------------------------------------------------------------")
print ("                                 gerando boxplot                              ")
print ("                             de sequencias com outliers                       ")
print ("                                       HA                                     ")
print ("------------------------------------------------------------------------------")

todos_codons_com_out1 = todos_codons_com_out[todos_codons_com_out["Segment"] == "HA"]
com_out_sorted = todos_codons_com_out1[sorted(todos_codons_com_out1.columns)]
plt.title("Box Plot dos Códons de Sequências do Vírus Influenza A para o gene HA com Outliers")
plt.xlabel("Códons")
plt.ylabel("Contagem")
plt.xticks(rotation=90)
sns.boxplot(data=com_out_sorted, color="lightgray", boxprops=dict(facecolor='gray'),showfliers=False)
plt.show()
print(com_out_sorted)
com_out_sorted.info()

print("------------------------------------------------------------------------------")
print("                                 gerando boxplot                              ")
print("                             de sequencias com outliers                       ")
print("                                       NA                                     ")
print("------------------------------------------------------------------------------")

todos_codons_com_out2 = todos_codons_com_out[todos_codons_com_out["Segment"] == "NA"]
com_out_sorted2 = todos_codons_com_out2[sorted(todos_codons_com_out2.columns)]
plt.title("Box Plot dos Códons de Sequências do Vírus Influenza A para o gene NA com Outliers")
plt.xlabel("Códons")
plt.ylabel("Contagem")
plt.xticks(rotation=90)
sns.boxplot(data=com_out_sorted2, color="lightgray", boxprops=dict(facecolor='gray'),showfliers=False)
plt.show()
print(com_out_sorted2)
com_out_sorted2.info()

print("------------------------------------------------------------------------------")
print("                              Exibindo todos os codons                        ")
print("                             de sequencias sem outliers                       ")
print("------------------------------------------------------------------------------")
df_filtrado2 = df1[df1["ID"].isin(df_sem_out["ID"])].reset_index(drop=True)
todos_codons_sem_out = pd.concat([df_filtrado2[["ID", "Segment"]], codon_df.reset_index(drop=True)], axis=1)
todos_codons_sem_out = todos_codons_sem_out.dropna(subset=["ID"])
print(todos_codons_sem_out)

print("------------------------------------------------------------------------------")
print("                                 gerando boxplot                              ")
print("                             de sequencias sem outliers                       ")
print("                                       HA                                     ")
print("------------------------------------------------------------------------------")
todos_codons_sem_out1 = todos_codons_sem_out[todos_codons_sem_out["Segment"] == "HA"]
sem_out_sorted = todos_codons_sem_out1[sorted(todos_codons_sem_out1.columns)]
plt.title("Box Plot dos Códons de Sequências do Vírus Influenza A para o gene HA sem Outliers")
plt.xlabel("Códons")
plt.ylabel("Contagem")
plt.xticks(rotation=90)
sns.boxplot(data=sem_out_sorted, color="lightgray", boxprops=dict(facecolor='gray'),showfliers=False)
plt.show()
print(sem_out_sorted)

print("------------------------------------------------------------------------------")
print("                                 gerando boxplot                              ")
print("                             de sequencias sem outliers                       ")
print("                                       NA                                     ")
print("------------------------------------------------------------------------------")
todos_codons_sem_out2 = todos_codons_sem_out[todos_codons_sem_out["Segment"] == "NA"]
sem_out_sorted2 = todos_codons_sem_out2[sorted(todos_codons_sem_out2.columns)]
plt.title("Box Plot dos Códons de Sequências do Vírus Influenza A para o gene NA sem Outliers")
plt.xlabel("Códons")
plt.ylabel("Contagem")
plt.xticks(rotation=90)
sns.boxplot(data=sem_out_sorted2, color="lightgray", boxprops=dict(facecolor='gray'),showfliers=False)
plt.show()
print(sem_out_sorted2)

print("------------------------------------------------------------------------------")
print("                                 gerando gráfico                              ")
print("                  comparando sequencias com e sem outliers                    ")
print("                                     de HA                                    ")
print("------------------------------------------------------------------------------")
com_out_sorted = com_out_sorted.drop(columns=["ID", "Segment"])
sem_out_sorted = sem_out_sorted.drop(columns=["ID", "Segment"])
comparacao = pd.DataFrame({
    "Não outlier": sem_out_sorted.mean(),
    "Outlier": com_out_sorted.mean()
})

comparacao.plot(kind="bar", color=["gray", "black"], alpha=0.8)
plt.xlabel("Códons")
plt.ylabel("Contagem Média")
plt.title("Comparação de Média da Quantidade de Códons de Sequências com e sem Outliers para o gene HA")
plt.xticks(rotation=90)
plt.legend()
plt.tight_layout()
plt.show()
# Gráfico de dispersão das médias#
codons = comparacao.index
x_pos = range(len(codons))

plt.figure(figsize=(12, 6))

# Plotando os pontos de Contagem Média para "Não outlier"
plt.scatter(x_pos, comparacao["Não outlier"], color="gray", label="Não outlier", marker="o", s=50)

# Plotando os pontos de Contagem Média para "Outlier"
plt.scatter(x_pos, comparacao["Outlier"], color="black", label="Outlier", marker="x", s=50)

plt.xlabel("Códons")
plt.ylabel("Contagem Média")
plt.title("Comparação de Média da Quantidade de Códons de Sequências com e sem Outliers para o gene HA (Dispersão de Médias)")
plt.xticks(x_pos, codons, rotation=90) # Define os rótulos do eixo X com os nomes dos códons
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6) # Adiciona linhas de grade para melhor leitura
plt.tight_layout()
plt.show()


print("------------------------------------------------------------------------------")
print("                                 gerando gráfico                              ")
print("                  comparando sequencias com e sem outliers                    ")
print("                                     de NA                                    ")
print("------------------------------------------------------------------------------")
com_out_sorted2 = com_out_sorted2.drop(columns=["ID", "Segment"])
sem_out_sorted2 = sem_out_sorted2.drop(columns=["ID", "Segment"])
comparacao2 = pd.DataFrame({
    "Não outlier": sem_out_sorted2.mean(),
    "Outlier": com_out_sorted2.mean()
})

comparacao2.plot(kind="bar", color=["gray", "black"], alpha=0.8)
plt.xlabel("Códons")
plt.ylabel("Contagem Média")
plt.title("Comparação de Média da Quantidade de Códons de Sequências com e sem Outliers para o gene NA")
plt.xticks(rotation=90)
plt.legend()
plt.tight_layout()
plt.show()
# Gráfico de dispersão das médias#
codons2 = comparacao2.index
x_pos2 = range(len(codons2))

plt.figure(figsize=(12, 6))

# Plotando os pontos de Contagem Média para "Não outlier"
plt.scatter(x_pos2, comparacao2["Não outlier"], color="gray", label="Não outlier", marker="o", s=50)

# Plotando os pontos de Contagem Média para "Outlier"
plt.scatter(x_pos2, comparacao2["Outlier"], color="black", label="Outlier", marker="x", s=50)

plt.xlabel("Códons")
plt.ylabel("Contagem Média")
plt.title("Comparação de Média da Quantidade de Códons de Sequências com e sem Outliers para o gene NA (Dispersão de Médias)")
plt.xticks(x_pos2, codons2, rotation=90)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout



# Separando os dataframes SEM outliers por segmento (HA e NA)
df_ha_sem_out = todos_codons_sem_out[todos_codons_sem_out['Segment'] == 'HA'].drop(columns=['ID', 'Segment'])
df_na_sem_out = todos_codons_sem_out[todos_codons_sem_out['Segment'] == 'NA'].drop(columns=['ID', 'Segment'])

# Calculando a mediana das frequências de códons
mediana_ha_sem_out = df_ha_sem_out.median(numeric_only=True)
mediana_na_sem_out = df_na_sem_out.median(numeric_only=True)

# Unindo as medianas em um único DataFrame e ordenando por HA
df_medianas_sem_out = pd.DataFrame({'HA': mediana_ha_sem_out, 'NA': mediana_na_sem_out})
df_medianas_sem_out = df_medianas_sem_out.sort_values(by='HA', ascending=True)

plt.figure(figsize=(12, 8))
sns.scatterplot(x='HA', y='NA', data=df_medianas_sem_out)
plt.plot([0, df_medianas_sem_out['HA'].max()], [0, df_medianas_sem_out['HA'].max()], 'r--', label='y=x (Correlação Perfeita)')

# Adicionando rótulos aos pontos
for i, txt in enumerate(df_medianas_sem_out.index):
    plt.annotate(txt, (df_medianas_sem_out['HA'].iloc[i], df_medianas_sem_out['NA'].iloc[i]), xytext=(-4, 4),textcoords='offset points',fontsize=8, alpha=0.7)

plt.title('Correlação do Uso de Códons (HA vs NA) - Sem Outliers')
plt.xlabel('Frequência Relativa Mediana (Gene HA)')
plt.ylabel('Frequência Relativa Mediana (Gene NA)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

# --- Gráfico de correlação COM Outliers ---

# Separando os dataframes COM outliers por segmento (HA e NA)
df_ha_com_out = todos_codons_com_out[todos_codons_com_out['Segment'] == 'HA'].drop(columns=['ID', 'Segment'])
df_na_com_out = todos_codons_com_out[todos_codons_com_out['Segment'] == 'NA'].drop(columns=['ID', 'Segment'])

# Calculando a mediana das frequências de códons
mediana_ha_com_out = df_ha_com_out.median(numeric_only=True)
mediana_na_com_out = df_na_com_out.median(numeric_only=True)

# Unindo as medianas em um único DataFrame e ordenando por HA
df_medianas_com_out = pd.DataFrame({'HA': mediana_ha_com_out, 'NA': mediana_na_com_out})
df_medianas_com_out = df_medianas_com_out.sort_values(by='HA', ascending=True)

plt.figure(figsize=(12, 8))
sns.scatterplot(x='HA', y='NA', data=df_medianas_com_out)
plt.plot([0, df_medianas_com_out['HA'].max()], [0, df_medianas_sem_out['HA'].max()], 'r--', label='y=x (Correlação Perfeita)')

# Adicionando rótulos aos pontos
for i, txt in enumerate(df_medianas_com_out.index):
    plt.annotate(txt, (df_medianas_com_out['HA'].iloc[i], df_medianas_com_out['NA'].iloc[i]),xytext=(-4, 4),textcoords='offset points', fontsize=8, alpha=0.7)

plt.title('Correlação do Uso de Códons (HA vs NA) - Com Outliers')
plt.xlabel('Frequência Relativa Mediana (Gene HA)')
plt.ylabel('Frequência Relativa Mediana (Gene NA)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()




df_comparacao_ha = pd.DataFrame({
    'HA Sem Outlier': mediana_ha_sem_out,
    'HA Com Outlier': mediana_ha_com_out
})
df_comparacao_ha = df_comparacao_ha.sort_values(by='HA Sem Outlier', ascending=True)

plt.figure(figsize=(12, 8))
sns.scatterplot(x='HA Sem Outlier', y='HA Com Outlier', data=df_comparacao_ha)

# Adiciona a linha de referência y=x
max_val_ha = max(df_comparacao_ha['HA Sem Outlier'].max(), df_comparacao_ha['HA Com Outlier'].max())
plt.plot([0, max_val_ha], [0, max_val_ha], 'r--', label='y=x (Nenhum Impacto do Outlier)')

# Adicionando rótulos aos pontos (nomes dos códons)
for i, txt in enumerate(df_comparacao_ha.index):
    plt.annotate(
        txt, 
        (df_comparacao_ha['HA Sem Outlier'].iloc[i], df_comparacao_ha['HA Com Outlier'].iloc[i]), 
        xytext=(-4, 4), 
        textcoords='offset points', 
        fontsize=8, 
        alpha=0.7
    )

plt.title('Impacto de Outliers na Mediana de Frequência de Códons (Gene HA)')
plt.xlabel('Frequência Relativa Mediana (HA Sem Outlier)')
plt.ylabel('Frequência Relativa Mediana (HA Com Outlier)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()



df_comparacao_na = pd.DataFrame({
    'NA Sem Outlier': mediana_na_sem_out,
    'NA Com Outlier': mediana_na_com_out
})
df_comparacao_na = df_comparacao_na.sort_values(by='NA Sem Outlier', ascending=True)

plt.figure(figsize=(12, 8))
sns.scatterplot(x='NA Sem Outlier', y='NA Com Outlier', data=df_comparacao_na)

# Adiciona a linha de referência y=x
max_val_na = max(df_comparacao_na['NA Sem Outlier'].max(), df_comparacao_na['NA Com Outlier'].max())
plt.plot([0, max_val_na], [0, max_val_na], 'r--', label='y=x (Nenhum Impacto do Outlier)')

# Adicionando rótulos aos pontos (nomes dos códons)
for i, txt in enumerate(df_comparacao_na.index):
    plt.annotate(
        txt, 
        (df_comparacao_na['NA Sem Outlier'].iloc[i], df_comparacao_na['NA Com Outlier'].iloc[i]), 
        xytext=(-4, 4), 
        textcoords='offset points', 
        fontsize=8, 
        alpha=0.7
    )

plt.title('Impacto de Outliers na Mediana de Frequência de Códons (Gene NA)')
plt.xlabel('Frequência Relativa Mediana (NA Sem Outlier)')
plt.ylabel('Frequência Relativa Mediana (NA Com Outlier)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()






print ("------------------------------------------------------------------------------")
print ("                                Gerando arquivo fasta                         ")
print ("                      Para  analisar sequências com outliers                  ")
print ("------------------------------------------------------------------------------")
sequencias_com_out_fasta = df1[df1['ID'].isin(df_com_out['ID'])]
output_path = os.path.join(download_path, "filtrado.fa")
with open(output_path, "w") as f: 
    for _, row in sequencias_com_out_fasta.iterrows():
        f.write(f">{row['ID']}\n{row['Sequence']}\n")
print(f"Arquivo FASTA salvo em: {output_path}")