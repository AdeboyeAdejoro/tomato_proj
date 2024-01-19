import pandas as pd
import os
from tqdm import tqdm

dirs = ["Lyc", "Pen"]
lyc_dataframes = {}
pen_dataframes = {}

for directory in dirs:
    bed_dir_path = f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/{directory}/bedfile_folder"
    for bed_file in tqdm(os.listdir(bed_dir_path)):
        print(f"starting loop for {bed_file}")
        bed_file_path = os.path.join(bed_dir_path, bed_file)
        columns = ["chromosome", "start", "end", "name", "score", "strand"]
        keyname = f"{bed_file}_dataframe"
        if directory == "Lyc":
            lyc_dataframes[keyname] = pd.read_csv(bed_file_path, sep="\t", header=None, names=columns)
        elif directory == "Pen":
            pen_dataframes[keyname] = pd.read_csv(bed_file_path, sep="\t", header=None, names=columns)


common_keys = set(lyc_dataframes.keys()).intersection(pen_dataframes.keys())
print("common keys obtained")

common_names_dict = {}
lyc_unique_dict = {}
pen_unique_dict = {}

for key in common_keys:
    print(f"processing for {key}")
    lyc_names = set(lyc_dataframes[key]['name'])
    pen_names = set(pen_dataframes[key]['name'])
    
    common_names = lyc_names.intersection(pen_names)
    common_names_dict[key] = common_names
    
    lyc_unique_dict[key] = lyc_names - common_names
    pen_unique_dict[key] = pen_names - common_names


if not os.path.exists(f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files"):
    os.mkdir(f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files")

written_files_folder = f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files"

print("Attempting to write files")

common_names_df = pd.DataFrame.from_dict(common_names_dict, orient="index")
common_names_df.to_csv(f"{written_files_folder}/common_names_df.csv", index_label="key")

lyc_unique_df = pd.DataFrame.from_dict(lyc_unique_dict, orient="index")
lyc_unique_df.to_csv(f"{written_files_folder}/lyc_unique_df.csv", index_label="key")

pen_unique_df = pd.DataFrame.from_dict(pen_unique_dict, orient="index")
pen_unique_df.to_csv(f"{written_files_folder}/pen_unique_df.csv", index_label="key")

