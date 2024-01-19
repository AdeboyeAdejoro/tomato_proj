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

union_names_dict = {}

for key in common_keys:
    print(f"processing for {key}")
    lyc_names = set(lyc_dataframes[key]['name'])
    pen_names = set(pen_dataframes[key]['name'])
    
    union_names = lyc_names.union(pen_names)
    union_names_dict[key] = union_names
    print("union_names_dictionary written")

if not os.path.exists(f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files_restructured"):
    os.mkdir(f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files_restructured")

# New folder for restructured files
written_files_folder_restructured = f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files_restructured"

#print("Converting union names dictionary to dataframe")
#union_names_df = pd.DataFrame.from_dict(union_names_dict, orient="index")


for key in tqdm(common_keys):
    print(f"Restructuring data for {key}")

    # Create a new DataFrame for each bed file
    current_df = pd.DataFrame()

    # Set the 'reads' column based on the union names
    current_df['reads'] = union_names_dict[key]

    # Create new columns for lycopersicum and pennelli, fill with 0 initially
    current_df[['lycopersicum', 'pennelli']] = 0

    # Mark reads that match lycopersicum with 1
    current_df['lycopersicum'] = current_df['reads'].isin(lyc_dataframes[key]['name']).astype(int)

    # Mark reads that match pennelli with 1
    current_df['pennelli'] = current_df['reads'].isin(pen_dataframes[key]['name']).astype(int)

    # Save the restructured dataframe to a CSV file
    current_df.to_csv(f"{written_files_folder_restructured}/{key}.csv", index=False)

print("Data restructuring completed.")

