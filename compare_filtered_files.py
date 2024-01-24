import os
import pandas as pd
from tqdm import tqdm

folder_path = "/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files_restructured"

if not os.path.exists("/nam-99/ablage/nam/adejoro_folder/ILSolanum/filtered"):
    os.mkdir("/nam-99/ablage/nam/adejoro_folder/ILSolanum/filtered")

filenames = []
empty_flags = []
not_empty_flags = []

for file in tqdm(os.listdir(folder_path)):
    print(f"looping over: {file}")
    file_path = os.path.join(folder_path, file)

    df_inner = pd.read_csv(file_path)
    df_inner['sum'] = df_inner['lycopersicum'] + df_inner['pennelli']

    filtered_df_inner = df_inner[(df_inner['sum'] != 2)]
    filtered_df_inner = df_inner[(df_inner['lycopersicum'] != 1)]

    filenames.append(file)
    empty_flags.append(int(filtered_df_inner.empty))
    not_empty_flags.append(int(not filtered_df_inner.empty))

df = pd.DataFrame({'filenames': filenames, 'empty': empty_flags, 'not empty': not_empty_flags})
df.to_csv("/nam-99/ablage/nam/adejoro_folder/ILSolanum/filtered/compare.csv", index=False)
print("Script finished")
