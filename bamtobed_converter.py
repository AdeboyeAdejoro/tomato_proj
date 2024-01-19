import pandas as pd
import os
from tqdm import tqdm
#pd.options.display.width = 0

dirs = ["Lyc", "Pen"]
#lyc_dataframes = {}
#pen_dataframes = {}

for directory in dirs:
    dir_path = f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/{directory}"
    if not os.path.exists(f"{dir_path}/bedfile_folder"):
        os.mkdir(f"{dir_path}/bedfile_folder")

    print(f"processing files in {dir_path}")

    for file in tqdm(os.listdir(dir_path)):
        print(f"loop for {file}")
        if os.path.isdir(os.path.join(dir_path, file)) or not file.endswith(".bam"):
            print(f"Skipping {file} as it is a folder or non bamfile")
            continue
        if file.endswith(".bam"):
            print(f"converting {file} to {file}.bed") 
            if not os.path.exists(f"{dir_path}/bedfile_folder/{file}.bed"):
                os.system(f"bedtools bamtobed -i {dir_path}/{file} > {dir_path}/bedfile_folder/{file}.bed")





