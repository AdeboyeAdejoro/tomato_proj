import os
import pandas as pd
from tqdm import tqdm

folder_path = "/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files_restructured"

if not os.path.exists("/nam-99/ablage/nam/adejoro_folder/ILSolanum/filtered"):
    os.mkdir("/nam-99/ablage/nam/adejoro_folder/ILSolanum/filtered")

selected_files = []

for file in tqdm(os.listdir(folder_path)):
    file_path = os.path.join(folder_path, file)
    
    print(f"reading {file}")
    df = pd.read_csv(file_path)
    
    df['sum_column'] = df['lycopersicum'] + df['pennelli']
    
    filtered_df = df[(df['sum_column'] !=2) & (df['pennelli'] != 1)]
    
    if not filtered_df.empty:
        selected_files.append(file)
        
result_df = pd.DataFrame({'Selected_Files' : selected_files})

result_df.to_csv("/nam-99/ablage/nam/adejoro_folder/ILSolanum/filtered/selected_files.csv", index=False)

print("File selection completed")
