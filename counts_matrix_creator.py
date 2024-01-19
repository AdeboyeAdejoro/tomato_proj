import os
import pandas as pd
from tqdm import tqdm

dirs = ['Lyc', 'Pen']

for dir in dirs:
    path_var = f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/{dir}/counts_folder"
    output_folder = f"/nam-99/ablage/nam/adejoro_folder/ILSolanum/{dir}/counts_matrix_folder"

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    final_df = pd.DataFrame()

    for file in tqdm(os.listdir(path_var)):

        file_path = os.path.join(path_var, file)

        # Skip empty files
        if os.path.getsize(file_path) == 0:
            print(f"Skipping empty file: {file}")
            continue

        # Read count file and process data
        read_count_file = pd.read_csv(file_path, sep="\t", comment="#")
        req_data = read_count_file.iloc[:, [0, -1]]
        transpose_req_data = req_data.T
        transpose_req_data.columns = transpose_req_data.iloc[0]
        transpose_req_data = transpose_req_data[1:]

        # Remove 'gene:' prefix from gene IDs
        transpose_req_data.columns = transpose_req_data.columns.str.replace('gene:', '')

        final_df = pd.concat([final_df, transpose_req_data], axis=0)

# Rename 'Geneid' column to 'SRR_accession'
    final_df.columns = ['SRR_accession' if col == 'Geneid' else col for col in final_df.columns]

# Write the final dataframe to a CSV file
    output_file_path = os.path.join(output_folder, 'counts_matrix_with_index.csv')
    final_df.to_csv(output_file_path, index=True)

# Display a message
print(f"Final dataframe written")
