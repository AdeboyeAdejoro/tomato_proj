import os
from tqdm import tqdm


bamfiles_dir_path = "/nam-99/ablage/nam/adejoro_folder/ILSolanum/Lyc"
bamfiles_folder = "/nam-99/ablage/nam/adejoro_folder/bamfiles_folder"

if not os.path.exists(bamfiles_folder):	
	os.mkdir(bamfiles_folder)

for file in tqdm(os.listdir(bamfiles_dir_path)):
    constructed_file_path = os.path.join(bamfiles_dir_path, file)
    if os.path.isdir(constructed_file_path) or not file.endswith(".bam"):
            print(f"skipping folder or non bamfile")
            continue
    if os.path.isfile(constructed_file_path) and file.endswith(".bam"):
            os.system(f"cp {constructed_file_path} {bamfiles_folder}")

