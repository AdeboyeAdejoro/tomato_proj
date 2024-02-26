import os
import pandas as pd
from tqdm import tqdm
source_path = "/nam-99/ablage/nam/adejoro_folder/filtered_feature_counts/filtered_counts_only_folder"
dest_path = "/nam-99/ablage/nam/adejoro_folder/filtered_counts_matrix"

if not os.path.exists(dest_path):
    os.mkdir(dest_path)

final_df = pd.DataFrame()

for file in tqdm(os.listdir(source_path)):
    print(f"loop for file: {file}") 
    constructed_path = os.path.join(source_path, file)
    read_count_file = pd.read_csv(constructed_path, sep="\t", comment="#")
    req_data = read_count_file.iloc[:, [0, -1]]
    transpose_req_data = req_data.T
    transpose_req_data.columns = transpose_req_data.iloc[0]
    transpose_req_data = transpose_req_data[1:]
    transpose_req_data.columns = transpose_req_data.columns.str.replace('gene', '')
    final_df = pd.concat([final_df, transpose_req_data], axis=0)

output_path = os.path.join(dest_path, "filtered_counts_matrix.csv")
final_df.to_csv(output_path, index=True)
print("complete")
