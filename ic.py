from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
import time
import os


download_path = r"D:\codigopython\ArqFasta"  

# Configurar opções do Chrome
chrome_options = Options()
chrome_options.add_experimental_option("prefs", {
    "download.default_directory": download_path,  
    "download.prompt_for_download": False,       
    "download.directory_upgrade": True,         
    "safebrowsing.enabled": True                
})

# Configurar o ChromeDriver
service = Service(ChromeDriverManager().install())
driver = webdriver.Chrome(service=service, options=chrome_options)

def is_download_complete(download_dir):
    # Listar os arquivos no diretório de download
    files = os.listdir(download_dir)
    for file in files:
        file_path = os.path.join(download_dir, file)
        # Verifica se o arquivo tem o sufixo .crdownload (significa que não foi totalmente baixado) (na hora parece não funcionar mas reduz muito os erros)
        if file.endswith('.crdownload'):
            return False
    return True

def get_file(bf):
    time.sleep(1)
    after_files = set(os.listdir(download_path))
    new_files = after_files - bf
    if new_files:
        return list(new_files)[0]
    return None


# Abrir a página do NCBI
url = "https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi#mainform"
driver.get(url)




# Aguarde até que o botão de nucleotídeo esteja presente e clique nele
try:
    nucleotide_button = WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.XPATH, '//input[@name="sequence" and @value="N"]'))
    )
    nucleotide_button.click()
    print("Botão de nucleotídeo clicado com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao clicar no botão de nucleotídeo: {e}")

# Selecione "Type A"
try:
    select_type = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, 'type'))
    )
    select = Select(select_type)
    select.select_by_visible_text('A')
    print("Tipo 'A' selecionado com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao selecionar o tipo 'A': {e}")

# Selecione "Avian" no dropdown de Host
try:
    select_host = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, 'host'))
    )
    select = Select(select_host)
    select.select_by_visible_text('Avian')
    print("Host 'Avian' selecionado com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao selecionar o host 'Avian': {e}")

# Selecione "HA" no dropdown de Segment
try:
    select_segment = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, 'segment'))
    )
    select = Select(select_segment)
    
    select.select_by_value('4')
    print("Segmento 'HA' selecionado com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao selecionar o segmento 'HA': {e}")

x=0
def aposseg():
    
# Clique no checkbox "Full-length only"
    try:
        full_length_checkbox = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '//input[@name="sonly"]'))
        )
        full_length_checkbox.click()
        print("Checkbox 'Full-length only' clicado com sucesso!")
    except Exception as e:
        print(f"Ocorreu um erro ao clicar no checkbox 'Full-length only': {e}")

    # Selecione o checkbox "reqseg_full"
    try:
        reqseg_full_checkbox = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '//input[@name="reqseg_full" and @value="complete"]'))
        )
        reqseg_full_checkbox.click()
        print("Checkbox 'reqseg_full' clicado com sucesso!")
    except Exception as e:
        print(f"Ocorreu um erro ao clicar no checkbox 'reqseg_full': {e}")

    # Selecione o checkbox "reqseg"
    try:
        reqseg_checkbox = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '//input[@name="reqseg" and @value="4"]'))
        )
        reqseg_checkbox.click()
        print("Checkbox 'reqseg' clicado com sucesso!")
    except Exception as e:
        print(f"Ocorreu um erro ao clicar no checkbox 'reqseg': {e}")

    # Selecione a opção "Exclude" no dropdown "vac_strain"
    try:
        select_vac_strain = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'vac_strain'))
        )
        select = Select(select_vac_strain)
        select.select_by_visible_text('Exclude')
        print("Opção 'Exclude' selecionada com sucesso no dropdown 'vac_strain'!")
    except Exception as e:
        print(f"Ocorreu um erro ao selecionar no dropdown 'vac_strain': {e}")

    # Selecione a opção "Exclude" no dropdown "subtype_mix"
    try:
        select_subtype_mix = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'subtype_mix'))
        )
        select = Select(select_subtype_mix)
        select.select_by_visible_text('Exclude')
        print("Opção 'Exclude' selecionada com sucesso no dropdown 'subtype_mix'!")
    except Exception as e:
        print(f"Ocorreu um erro ao selecionar no dropdown 'subtype_mix': {e}")

    # Loop para H de H1 até H18
    for h in range(1, 19):
        WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'subtype_n'))
        )
        subtype_h = str(h)
        
        # Loop para N de N1 até N11
        for n in range(1, 12):
            subtype_n = str(n)
            bf = set(os.listdir(download_path))
            # Selecionar H
            try:
                select_h = Select(driver.find_element(By.NAME, 'subtype_h'))
                select_h.deselect_all()
                select_h.select_by_visible_text(subtype_h)
                print(f"Subtipo H '{subtype_h}' selecionado com sucesso!")

                # Selecionar N
                select_n = Select(driver.find_element(By.NAME, 'subtype_n'))
                select_n.deselect_all()
                select_n.select_by_visible_text(subtype_n)
                print(f"Subtipo N '{subtype_n}' selecionado com sucesso!")

                # Clique no botão 'Show results'
                show_results_button = WebDriverWait(driver, 10).until(
                    EC.element_to_be_clickable((By.XPATH, '//span[@class="default combined"]//input[@name="cmd_get_query"]'))
                )
                show_results_button.click()
                print(f"Botão 'Show results' clicado para H{subtype_h}N{subtype_n}!")

                # Aguarde um pouco entre as iterações
                time.sleep(1)
                windows = driver.window_handles  # Captura todas as janelas abertas
                if len (windows)>1:
                    driver.switch_to.window(windows[-1])  # Muda para a última janela (a nova)
                # Agora, clique no botão de "Download"
                # Aguarde até que o botão de "Download" esteja visível e habilitado
                try:
                    download_button = WebDriverWait(driver, 30).until(
                        EC.visibility_of_element_located((By.NAME, 'download-button'))
                    )
                    download_button.click()
                    
                    WebDriverWait(driver, 120, 1).until(lambda driver: is_download_complete(download_path))
                    time.sleep(6)
                    print("Botão de 'Download' clicado com sucesso!")
                    
                    new_file = get_file(bf)
                    if new_file:
                        new_file_path = os.path.join(download_path, new_file)
                        suffix = "_HA_" if x == 0 else "_NA_" if x == 1 else ""
                        renamed_file_path = os.path.join(download_path, f"H{subtype_h}N{subtype_n}(avian){suffix}.fa")
                        os.rename(new_file_path, renamed_file_path)
                        print(f"Arquivo renomeado")
                except Exception as e:
                    print(f"Erro ao clicar no botão de 'Download': {e}")


                # Aguarde um pouco após o clique
                time.sleep(1)
                if len (windows)>2:
                    
                    driver.close()  # Volta para a página onde as seleções foram feitas
                    windows = driver.window_handles  # Captura todas as janelas abertas
                    driver.switch_to.window(windows[-1])  # Muda para a última janela (a nova)
                    # Aguarde um pouco após voltar
                    time.sleep(1)
                
                else:
                    driver.back()
                
            except Exception as e:
                print(f"Ocorreu um erro ao processar H{subtype_h}N{subtype_n}: {e}")

    # Após o loop de "Avian", desmarque a seleção de "Avian" e selecione "Human"
    try:
        # Desmarque a opção "Avian"
        select_host = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'host'))
        )
        select = Select(select_host)
        select.deselect_all()  # Desmarcar todas as opções

        # Selecione "Human" no dropdown de Host
        select.select_by_visible_text('Human')  # Seleciona a opção "Human"
        print("Host 'Human' selecionado com sucesso!")

        # Repetir o loop de H e N para "Human"
        for h in range(1, 19):
            subtype_h = str(h)
            
            for n in range(1, 12):
                subtype_n = str(n)
                bf = set(os.listdir(download_path))
                # Selecionar H
                try:
                    select_h = Select(driver.find_element(By.NAME, 'subtype_h'))
                    select_h.deselect_all()
                    select_h.select_by_visible_text(subtype_h)
                    print(f"Subtipo H '{subtype_h}' selecionado com sucesso!")

                    # Selecionar N
                    select_n = Select(driver.find_element(By.NAME, 'subtype_n'))
                    select_n.deselect_all()
                    select_n.select_by_visible_text(subtype_n)
                    print(f"Subtipo N '{subtype_n}' selecionado com sucesso!")

                    # Clique no botão 'Show results'
                    show_results_button = WebDriverWait(driver, 10).until(
                        EC.element_to_be_clickable((By.XPATH, '//span[@class="default combined"]//input[@name="cmd_get_query"]'))
                    )
                    show_results_button.click()
                    print(f"Botão 'Show results' clicado para H{subtype_h}N{subtype_n}!")

                    # Aguarde um pouco entre as iterações
                    time.sleep(1)
                    windows = driver.window_handles  # Captura todas as janelas abertas
                    if len (windows)>2:
                        driver.switch_to.window(windows[-1])  # Muda para a última janela (a nova)
                    # Agora, clique no botão de "Download"
                    # Aguarde até que o botão de "Download" esteja visível e habilitado
                    try:
                        download_button = WebDriverWait(driver, 30).until(
                            EC.visibility_of_element_located((By.NAME, 'download-button'))
                        )
                        download_button.click()
                        WebDriverWait(driver, 120, 1).until(lambda driver: is_download_complete(download_path))
                        time.sleep(6)
                        
                        print("Botão de 'Download' clicado com sucesso!")
                        
                        new_file = get_file(bf)
                        if new_file:
                            new_file_path = os.path.join(download_path, new_file)
                            suffix = "_HA_" if x == 0 else "_NA_" if x == 1 else ""
                            renamed_file_path = os.path.join(download_path, f"H{subtype_h}N{subtype_n}(human){suffix}.fa")
                            os.rename(new_file_path, renamed_file_path)

                            #teste
                            
                            
                    except Exception as e:
                        print(f"Erro ao clicar no botão de 'Download': {e}")


                    # Aguarde um pouco após o clique
                    time.sleep(1)
                    if len (windows)>2:
                        
                        driver.close()  # Volta para a página onde as seleções foram feitas
                        windows = driver.window_handles  # Captura todas as janelas abertas
                        driver.switch_to.window(windows[-1])  # Muda para a última janela (a nova)
                        # Aguarde um pouco após voltar
                        time.sleep(1)
                    
                    else:
                        driver.back()

                except Exception as e:
                    print(f"Ocorreu um erro ao processar H{subtype_h}N{subtype_n}: {e}")

    except Exception as e:
        print(f"Ocorreu um erro ao selecionar o host 'Human': {e}")

    # Após o loop de "Human", desmarque a seleção de "Human" e selecione "Swine"
    try:
        # Desmarque a opção "Human"
        select_host = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'host'))
        )
        select = Select(select_host)
        select.deselect_all()  # Desmarcar todas as opções

        # Selecione "Swine" no dropdown de Host
        select.select_by_visible_text('Swine')  # Seleciona a opção "Swine"
        print("Host 'Swine' selecionado com sucesso!")

        # Repetir o loop de H e N para "Swine"
        for h in range(1, 19):#h19n12
            subtype_h = str(h)
            
            for n in range(1, 12):
                subtype_n = str(n)
                bf = set(os.listdir(download_path))
                # Selecionar H
                try:
                    select_h = Select(driver.find_element(By.NAME, 'subtype_h'))
                    select_h.deselect_all()
                    select_h.select_by_visible_text(subtype_h)
                    print(f"Subtipo H '{subtype_h}' selecionado com sucesso!")

                    # Selecionar N
                    select_n = Select(driver.find_element(By.NAME, 'subtype_n'))
                    select_n.deselect_all()
                    select_n.select_by_visible_text(subtype_n)
                    print(f"Subtipo N '{subtype_n}' selecionado com sucesso!")

                    # Clique no botão "Show results"
                    show_results_button = WebDriverWait(driver, 10).until(
                        EC.element_to_be_clickable((By.XPATH, '//span[@class="default combined"]//input[@name="cmd_get_query"]'))
                    )
                    show_results_button.click()
                    print(f"Botão 'Show results' clicado para H{subtype_h}N{subtype_n}!")

                    # Aguarde um pouco entre as iterações
                    time.sleep(1)
                    windows = driver.window_handles  # Captura todas as janelas abertas
                    if len (windows)>2:
                        driver.switch_to.window(windows[-1])  # Muda para a última janela (a nova)
                    # Agora, clique no botão de "Download"
                    # Aguarde até que o botão de "Download" esteja visível e habilitado
                    try:
                        download_button = WebDriverWait(driver, 30).until(
                            EC.visibility_of_element_located((By.NAME, 'download-button'))
                        )
                        download_button.click()
                        WebDriverWait(driver, 120, 1).until(lambda driver: is_download_complete(download_path))
                        time.sleep(3)
                        print("Botão de 'Download' clicado com sucesso!")
                        
                        new_file = get_file(bf)
                        if new_file:
                            new_file_path = os.path.join(download_path, new_file)
                            suffix = "_HA_" if x == 0 else "_NA_" if x == 1 else ""
                            renamed_file_path = os.path.join(download_path, f"H{subtype_h}N{subtype_n}(swine){suffix}.fa")
                            os.rename(new_file_path, renamed_file_path)
                            
                            print(f"Arquivo renomeado")
                    except Exception as e:
                        print(f"Erro ao clicar no botão de 'Download': {e}")


                    # Aguarde um pouco após o clique
                    time.sleep(1)
                    if len (windows)>2:
                        
                        driver.close()  # Volta para a página onde as seleções foram feitas
                        windows = driver.window_handles  # Captura todas as janelas abertas
                        driver.switch_to.window(windows[-1])  # Muda para a última janela (a nova)
                        # Aguarde um pouco após voltar
                        time.sleep(1)
                    
                    else:
                        driver.back()

                except Exception as e:
                    print(f"Ocorreu um erro ao processar H{subtype_h}N{subtype_n}: {e}")

    except Exception as e:
        print(f"Ocorreu um erro ao selecionar o host 'Swine': {e}")
aposseg()
try:
    select_segment = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, 'segment'))
    )
    select = Select(select_segment)
    select.select_by_value('6')
    select.deselect_by_value('4')
    print("Segmento 'NA' selecionado com sucesso!") #esta marcando junto ha
    x =x +1
except Exception as e:
    print(f"Ocorreu um erro ao selecionar o segmento 'NA': {e}")


aposseg()
# Fechar o navegador
driver.quit()