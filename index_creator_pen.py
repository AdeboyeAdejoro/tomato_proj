import os
from tqdm import tqdm

input_path = "/nam-99/ablage/nam/adejoro_folder/bamfiles_folder_pen"

for file in tqdm(os.listdir(input_path)):
	print(f"Loop for file: {file}")
	constructed_path = os.path.join(input_path, file)
	os.system(f"samtools index {constructed_path}")
