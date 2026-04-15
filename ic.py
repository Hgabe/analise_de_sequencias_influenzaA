from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
import time
import os


download_path = r"D:\codigopython\PROGRESSDOCODIGO\ArqFasta"  


chrome_options = Options()
chrome_options.add_experimental_option("prefs", {
    "download.default_directory": download_path,  
    "download.prompt_for_download": False,       
    "download.directory_upgrade": True,         
    "safebrowsing.enabled": True                
})


service = Service(ChromeDriverManager().install())
driver = webdriver.Chrome(service=service, options=chrome_options)

def is_download_complete(download_dir):
 
    files = os.listdir(download_dir)
    for file in files:
        file_path = os.path.join(download_dir, file)
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


url = "https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi#mainform"
driver.get(url)





try:
    nucleotide_button = WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.XPATH, '//input[@name="sequence" and @value="N"]'))
    )
    nucleotide_button.click()
    print("Botão de nucleotídeo clicado com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao clicar no botão de nucleotídeo: {e}")


try:
    select_type = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, 'type'))
    )
    select = Select(select_type)
    select.select_by_visible_text('A')
    print("Tipo 'A' selecionado com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao selecionar o tipo 'A': {e}")


try:
    select_host = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.NAME, 'host'))
    )
    select = Select(select_host)
    select.deselect_all()
    select.select_by_visible_text('Avian')
    print("Host 'Avian' selecionado com sucesso!")
except Exception as e:
    print(f"Ocorreu um erro ao selecionar o host 'Avian': {e}")


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
    

    try:
        full_length_checkbox = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '//input[@name="sonly"]'))
        )
        full_length_checkbox.click()
        print("Checkbox 'Full-length only' clicado com sucesso!")
    except Exception as e:
        print(f"Ocorreu um erro ao clicar no checkbox 'Full-length only': {e}")

 
    try:
        reqseg_full_checkbox = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '//input[@name="reqseg_full" and @value="complete"]'))
        )
        reqseg_full_checkbox.click()
        print("Checkbox 'reqseg_full' clicado com sucesso!")
    except Exception as e:
        print(f"Ocorreu um erro ao clicar no checkbox 'reqseg_full': {e}")

    
    try:
        reqseg_checkbox = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '//input[@name="reqseg" and @value="4"]'))
        )
        reqseg_checkbox.click()
        print("Checkbox 'reqseg' clicado com sucesso!")
    except Exception as e:
        print(f"Ocorreu um erro ao clicar no checkbox 'reqseg': {e}")

   
    try:
        select_vac_strain = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'vac_strain'))
        )
        select = Select(select_vac_strain)
        select.select_by_visible_text('Exclude')
        print("Opção 'Exclude' selecionada com sucesso no dropdown 'vac_strain'!")
    except Exception as e:
        print(f"Ocorreu um erro ao selecionar no dropdown 'vac_strain': {e}")

 
    try:
        select_subtype_mix = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'subtype_mix'))
        )
        select = Select(select_subtype_mix)
        select.select_by_visible_text('Exclude')
        print("Opção 'Exclude' selecionada com sucesso no dropdown 'subtype_mix'!")
    except Exception as e:
        print(f"Ocorreu um erro ao selecionar no dropdown 'subtype_mix': {e}")

    try:
        select_host = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'host'))
        )
        select = Select(select_host)
        select.deselect_all() 
        select.select_by_visible_text('Avian') 
        print("Host 'Avian' REINICIADO com sucesso para o primeiro loop!")
    except Exception as e:
        print(f"Ocorreu um erro ao selecionar o host 'Avian': {e}")
    for h in range(1, 19):
        WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'subtype_n'))
        )
        subtype_h = str(h)
        
   
        for n in range(1, 12):
            subtype_n = str(n)
            bf = set(os.listdir(download_path))
 
            try:
                select_h = Select(driver.find_element(By.NAME, 'subtype_h'))
                select_h.deselect_all()
                select_h.select_by_visible_text(subtype_h)
                print(f"Subtipo H '{subtype_h}' selecionado com sucesso!")

          
                select_n = Select(driver.find_element(By.NAME, 'subtype_n'))
                select_n.deselect_all()
                select_n.select_by_visible_text(subtype_n)
                print(f"Subtipo N '{subtype_n}' selecionado com sucesso!")

               
                show_results_button = WebDriverWait(driver, 10).until(
                    EC.element_to_be_clickable((By.XPATH, '//span[@class="default combined"]//input[@name="cmd_get_query"]'))
                )
                show_results_button.click()
                print(f"Botão 'Show results' clicado para H{subtype_h}N{subtype_n}!")

               
                time.sleep(1)
                windows = driver.window_handles  
                if len (windows)>1:
                    driver.switch_to.window(windows[-1]) 
                
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

                time.sleep(1)
                if len (windows)>2:
                    
                    driver.close() 
                    windows = driver.window_handles 
                    driver.switch_to.window(windows[-1])  
              
                    time.sleep(1)
                
                else:
                    driver.back()
                
            except Exception as e:
                print(f"Ocorreu um erro ao processar H{subtype_h}N{subtype_n}: {e}")

    
    try:
        
        select_host = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'host'))
        )
        select = Select(select_host)
        select.deselect_all() 

        
        select.select_by_visible_text('Human') 
        print("Host 'Human' selecionado com sucesso!")

        
        for h in range(1, 19):
            subtype_h = str(h)
            
            for n in range(1, 12):
                subtype_n = str(n)
                bf = set(os.listdir(download_path))
         
                try:
                    select_h = Select(driver.find_element(By.NAME, 'subtype_h'))
                    select_h.deselect_all()
                    select_h.select_by_visible_text(subtype_h)
                    print(f"Subtipo H '{subtype_h}' selecionado com sucesso!")

                   
                    select_n = Select(driver.find_element(By.NAME, 'subtype_n'))
                    select_n.deselect_all()
                    select_n.select_by_visible_text(subtype_n)
                    print(f"Subtipo N '{subtype_n}' selecionado com sucesso!")

                    
                    show_results_button = WebDriverWait(driver, 10).until(
                        EC.element_to_be_clickable((By.XPATH, '//span[@class="default combined"]//input[@name="cmd_get_query"]'))
                    )
                    show_results_button.click()
                    print(f"Botão 'Show results' clicado para H{subtype_h}N{subtype_n}!")

             
                    time.sleep(1)
                    windows = driver.window_handles 
                    if len (windows)>2:
                        driver.switch_to.window(windows[-1])
                    
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

                       
                            
                            
                    except Exception as e:
                        print(f"Erro ao clicar no botão de 'Download': {e}")


                    
                    time.sleep(1)
                    if len (windows)>2:
                        
                        driver.close()  
                        windows = driver.window_handles  
                        driver.switch_to.window(windows[-1])
                        
                        time.sleep(1)
                    
                    else:
                        driver.back()

                except Exception as e:
                    print(f"Ocorreu um erro ao processar H{subtype_h}N{subtype_n}: {e}")

    except Exception as e:
        print(f"Ocorreu um erro ao selecionar o host 'Human': {e}")

    
    try:
        
        select_host = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.NAME, 'host'))
        )
        select = Select(select_host)
        select.deselect_all()  

        
        select.select_by_visible_text('Swine') 
        print("Host 'Swine' selecionado com sucesso!")

        for h in range(1, 19):
            subtype_h = str(h)
            
            for n in range(1, 12):
                subtype_n = str(n)
                bf = set(os.listdir(download_path))
           
                try:
                    select_h = Select(driver.find_element(By.NAME, 'subtype_h'))
                    select_h.deselect_all()
                    select_h.select_by_visible_text(subtype_h)
                    print(f"Subtipo H '{subtype_h}' selecionado com sucesso!")

                    
                    select_n = Select(driver.find_element(By.NAME, 'subtype_n'))
                    select_n.deselect_all()
                    select_n.select_by_visible_text(subtype_n)
                    print(f"Subtipo N '{subtype_n}' selecionado com sucesso!")

                    
                    show_results_button = WebDriverWait(driver, 10).until(
                        EC.element_to_be_clickable((By.XPATH, '//span[@class="default combined"]//input[@name="cmd_get_query"]'))
                    )
                    show_results_button.click()
                    print(f"Botão 'Show results' clicado para H{subtype_h}N{subtype_n}!")

                
                    time.sleep(1)
                    windows = driver.window_handles 
                    if len (windows)>2:
                        driver.switch_to.window(windows[-1]) 
                    
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


                    
                    time.sleep(1)
                    if len (windows)>2:
                        
                        driver.close()  
                        windows = driver.window_handles  
                        driver.switch_to.window(windows[-1])  
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
    print("Segmento 'NA' selecionado com sucesso!") 
    x =x +1
except Exception as e:
    print(f"Ocorreu um erro ao selecionar o segmento 'NA': {e}")


aposseg()

driver.quit()