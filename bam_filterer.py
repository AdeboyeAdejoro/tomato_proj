import os
from tqdm import tqdm
import pandas as pd
import pysam
from glob import glob

bamfile_dir_path = "/nam-99/ablage/nam/adejoro_folder/bamfiles_folder"
corresponding_table_dir_path = "/nam-99/ablage/nam/adejoro_folder/ILSolanum/written_files_restructured"
output_dir_path = "/nam-99/ablage/nam/adejoro_folder/filtered_bam_folder"

if not os.path.exists(output_dir_path):
    os.mkdir(output_dir_path)

bam_files = glob(os.path.join(bamfile_dir_path, "*.bam"))

for bam_file in tqdm(bam_files):
    if os.path.getsize(bam_file) == 0:
        print(f"skipping empty bam file: {bam_file}")
        continue
    print(f"current bam_file is: {bam_file}")
    csv_file = os.path.join(corresponding_table_dir_path, f"{os.path.basename(bam_file)}.bed_dataframe.csv")
    mapping_df = pd.read_csv(csv_file)
    mapping_table = dict(zip(mapping_df['reads'], zip(mapping_df['lycopersicum'], mapping_df['pennelli'])))
    interpolated_name = os.path.basename(bam_file)
    print(f"interpolated_name is: {interpolated_name}")
    second_level_interpolated_name = interpolated_name.replace(".bam", "_filtered.bam")
    print(f"second_level_interpolated name is {second_level_interpolated_name}")
    output_bamfile_path = os.path.join(output_dir_path, second_level_interpolated_name)
    print(f"output bamfile path is: {output_bamfile_path}")

    if not os.path.exists(output_bamfile_path):
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        output_bamfile = pysam.AlignmentFile(output_bamfile_path, "wb", template=bamfile)

        for read in bamfile.fetch():
            read_id = read.query_name

            if read_id in mapping_table and mapping_table[read_id] == (1, 0):
                output_bamfile.write(read)

        bamfile.close()
        output_bamfile.close()
    
