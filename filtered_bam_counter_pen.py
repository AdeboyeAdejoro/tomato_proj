from tqdm import tqdm
import os

source_folder = "/nam-99/ablage/nam/adejoro_folder/filtered_bam_folder_for_pen"
output_folder = "/nam-99/ablage/nam/adejoro_folder/filtered_feature_counts_pen"

#source_folder = "/nam-99/ablage/nam/adejoro_folder/round_two/filtered_bam_folder"
#output_folder = "/nam-99/ablage/nam/adejoro_folder/round_two/filtered_feature_counts"

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


for file in tqdm(os.listdir(source_folder)):
    print(f"Loop for file: {file}")

    constructed_path = os.path.join(source_folder, file)
    intermediate_name = os.path.basename(constructed_path)
    second_level_intermediate_name = intermediate_name.replace(".bam", "_counts_pen.txt")
    constructed_output_file_path = os.path.join(output_folder, second_level_intermediate_name)
    if not os.path.exists(constructed_output_file_path):
    	os.system(f"featureCounts -t mRNA -g Parent -a spenn_v2.0_gene_models_annot.gff -o {constructed_output_file_path} {constructed_path}")
