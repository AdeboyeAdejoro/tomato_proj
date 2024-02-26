from tqdm import tqdm
import os

source_folder = "/nam-99/ablage/nam/adejoro_folder/filtered_bam_folder"
output_folder = "/nam-99/ablage/nam/adejoro_folder/filtered_feature_counts"

#source_folder = "/nam-99/ablage/nam/adejoro_folder/round_two/filtered_bam_folder"
#output_folder = "/nam-99/ablage/nam/adejoro_folder/round_two/filtered_feature_counts"

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


for file in tqdm(os.listdir(source_folder)):
    print(f"Loop for file: {file}")

    constructed_path = os.path.join(source_folder, file)
    intermediate_name = os.path.basename(constructed_path)
    second_level_intermediate_name = intermediate_name.replace(".bam", "_counts.txt")
    constructed_output_file_path = os.path.join(output_folder, second_level_intermediate_name)
    if not os.path.exists(constructed_output_file_path):
    	os.system(f"featureCounts -t mRNA -g Parent -a ITAG4.0_gene_models.gff -o {constructed_output_file_path} {constructed_path}")
