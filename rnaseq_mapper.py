import pandas as pd
from tqdm import tqdm
import os

root_dir, species_names = ['Lyc', 'Pen'], ['Slyc', 'Spenn']
for rdir, species in zip(root_dir, species_names):
    if not os.path.exists(f'ILSolanum/{rdir}/hisat2_index/{species}.8.ht2'):
        chroms = ','.join([f"ILSolanum/{rdir}/CHROM/{x}" for x in os.listdir(f'ILSolanum/{rdir}/CHROM/')])
        os.system(f"hisat2-build -p 7 {chroms} ILSolanum/{rdir}/hisat2_index/{species}")

sra_meta_df = pd.read_csv('SRR_accessions.txt', header=None)
print(sra_meta_df.head())
for accession_id in tqdm(sra_meta_df[0]):
    if not os.path.exists(f'ILSolanum/Pen/{accession_id}_counts.txt') and not os.path.exists(f'ILSolanum/Lyc/{accession_id}_counts.txt'):
        os.system(f"/home/adejoro/anaconda3/envs/tomatoEnv/bin/fasterq-dump {accession_id} -p -e 30")
        # Trim with Sickle
        os.system(f"trimmomatic SE -phred33 {accession_id}.fastq {accession_id}_trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> trim_{accession_id}.log")
        os.system(f"fastqc {accession_id}_trim.fastq")

        # Align to Pennellii and Lycopersicum
        os.system(f"hisat2 -p 20 -x ILSolanum/Lyc/hisat2_index/Slyc -U {accession_id}_trim.fastq -S ILSolanum/Lyc/{accession_id}.sam")
        os.system(f"hisat2 -p 20 -x ILSolanum/Pen/hisat2_index/Spenn -U {accession_id}_trim.fastq -S ILSolanum/Pen/{accession_id}.sam")
        # Convert Sam to Bam
        # Lyc
        os.system(f"samtools view -b ILSolanum/Lyc/{accession_id}.sam > ILSolanum/Lyc/{accession_id}.unsorted.bam")
        os.system(f"samtools sort ILSolanum/Lyc/{accession_id}.unsorted.bam > ILSolanum/Lyc/{accession_id}.bam")
        os.system(f"featureCounts -t mRNA -g Parent -a ITAG4.0_gene_models.gff -o ILSolanum/Lyc/{accession_id}_counts.txt ILSolanum/Lyc/{accession_id}.bam")
        os.system("rm -rf ILSolanum/Lyc/*unsorted.bam")
        os.system("rm -rf ILSolanum/Lyc/*.sam")
        # Pen
        os.system(f"samtools view -b ILSolanum/Pen/{accession_id}.sam > ILSolanum/Pen/{accession_id}.unsorted.bam")
        os.system(f"samtools sort ILSolanum/Pen/{accession_id}.unsorted.bam > ILSolanum/Pen/{accession_id}.bam")
        os.system(f"featureCounts -t mRNA -g Parent -a spenn_v2.0_gene_models_annot.gff -o ILSolanum/Pen/{accession_id}_counts.txt ILSolanum/Pen/{accession_id}.bam")
        os.system("rm -rf ILSolanum/Pen/*unsorted.bam")
        os.system("rm -rf ILSolanum/Pen/*.sam")

        # clean
        os.system(f"cp *.zip ILSolanum/Lyc/")
        os.system(f"cp *.zip ILSolanum/Pen/")
        os.system(f"rm -rf *.zip")
        os.system(f"rm -rf *.fastq")
        os.system(f"rm -rf *.bam")
        os.system(f"rm -rf *.html")
