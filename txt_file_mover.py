import os
from tqdm import tqdm

dirs = ['Lyc', 'Pen']

for dir in dirs:
    source_path = f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/{dir}" 
    destination_path = f"{source_path}/counts_folder"
    
    if not os.path.exists(destination_path):
        os.mkdir(destination_path)
        
    for file in tqdm(os.listdir(source_path)):
        file_path = os.path.join(source_path, file)
        if os.path.isdir(file_path) or not file.endswith(".txt"):
            print("Skipping folder or non-count file")
            continue
        destination_file_path = os.path.join(destination_path, file)
        if file.endswith(".txt") and not os.path.exists(destination_file_path):
            print("Copying count file")
            os.system(f"cp {file_path} {destination_path}")
        

