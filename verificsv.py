import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

download_path = r"D:\codigopython\ArqFasta"

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
csv_path = os.path.join(download_path, "verif.csv")
df.to_csv(csv_path, index=False)
print(f"Arquivo CSV salvo em: {csv_path}")



print(df)
df.info()



