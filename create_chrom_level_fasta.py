from Bio import SeqIO
import os

if not os.path.exists('Spenn.fasta'):
    os.system("wget https://solgenomics.net/ftp/genomes/Solanum_pennellii/Spenn.fasta")

if not os.path.exists('S_lycopersicum_chromosomes.4.00.fa'):
    os.system("wget https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa")

if not os.path.exists("spenn_v2.0_gene_models_annot.gff"):
    os.system("wget https://solgenomics.net/ftp/genomes/Solanum_pennellii/spenn_v2.0_gene_models_annot.gff")

if not os.path.exists("ITAG4.0_gene_models.gff"):
    os.system("wget https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_gene_models.gff")

if not os.path.exists('ILSolanum'):
    os.mkdir('ILSolanum')

if not os.path.exists('ILSolanum/Lyc'):
    os.mkdir('ILSolanum/Lyc')
    os.mkdir('ILSolanum/Lyc/CHROM')
    os.mkdir('ILSolanum/Lyc/hisat2_index')

if not os.path.exists('ILSolanum/Pen'):
    os.mkdir('ILSolanum/Pen')
    os.mkdir('ILSolanum/Pen/CHROM')
    os.mkdir('ILSolanum/Pen/hisat2_index')


for species, filename, root_dir in zip(['Spenn', 'SLyc'], ['Spenn.fasta', 'S_lycopersicum_chromosomes.4.00.fa'], ['Pen', 'Lyc']):
    for rec in SeqIO.parse(filename, format='fasta'):
        if not os.path.exists(f'ILSolanum/{root_dir}/{rec.id}.fa'):
            SeqIO.write(SeqIO.SeqRecord(seq=rec.seq, id=rec.id), handle=f'ILSolanum/{root_dir}/CHROM/{rec.id}.fa', format='fasta')