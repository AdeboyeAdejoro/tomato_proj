import pandas as pd
import numpy as np
import os
from utils import FastaSequenceLoader, ConvNetwork
from tensorflow.keras.backend import clear_session
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
mapped_read_counts = ['pure_lyc_counts.csv']
genomes = ['S_lycopersicum_chromosomes.4.00.fa']
gene_models = ['ITAG4.0_gene_models.gff']
num_chromosomes = [['SL4.0ch01', 'SL4.0ch02', 'SL4.0ch03', 'SL4.0ch04', 'SL4.0ch05',
                  'SL4.0ch06', 'SL4.0ch07', 'SL4.0ch08', 'SL4.0ch09', 'SL4.0ch10', 'SL4.0ch11', 'SL4.0ch12']]

if not os.path.isdir('results'):
    os.mkdir('results')
if not os.path.isdir('saved_models'):
    os.mkdir('saved_models')

for m_reads, gene_model, genome, num_chr in zip(mapped_read_counts, gene_models, genomes, num_chromosomes):
    if not os.path.exists(f"results/{m_reads.split('_')[0]}_result.csv"):
        final_training_output = []
        counts = pd.read_csv(f"counts/{m_reads}", index_col=0)
        true_targets = []
        
        for log_count in counts['median'].values:
            """if log_count <=2:
                true_targets.append(0)
            else:
                true_targets.append(1)"""
            if log_count <= np.percentile(counts['median'], 25):
                true_targets.append(0)
            elif log_count >= np.percentile(counts['median'], 75):
                true_targets.append(1)
            else: 
                true_targets.append(2)
                
        counts['true_target'] = true_targets
        print(counts.head(20))
        
        for val_chromosome in num_chr:
            fastaloader = FastaSequenceLoader(f'genomes/{genome}', f'gene_models/{gene_model}', val_chromosome)
            enc_train, enc_val, train_ids, val_ids = fastaloader.extract_seq()
            print(len(enc_train), len(enc_val))
            print(train_ids[:10])
            print('---------------------------------------------------------------------\n')
            print(f"Plant: {m_reads.split('_')[0]} Case: promoter_terminator")
            print('-----------------------------------------------------------------------')
            convnet = ConvNetwork(enc_train, enc_val, train_ids, val_ids, val_chromosome, counts, m_reads.split('_')[0], 'promoter_terminator', tissue='fruit skin')
            clear_session()
            output = convnet.train_network()
            final_training_output.append(output)

            #Train models with shuffled sequences
            print('----------------------------------------------------------\n')
            print(f"Plant: {m_reads.split('_')[0]} Case: si-nucleotide_shuffle")
            print('-------------------------------------------------------------')
            shuffle_enc_train = []
            for train_seq in enc_train.copy():
                np.random.shuffle(train_seq)
                shuffle_enc_train.append(train_seq)

            shuffle_convnet = ConvNetwork(shuffle_enc_train, enc_val, train_ids, val_ids, val_chromosome, counts, m_reads.split('_')[0], 'si-nucleotide_shuffle', tissue="fruit skin")
            clear_session()
            shuffle_output = shuffle_convnet.train_network()
            print("\n before appending shuffle output")
            final_training_output.append(shuffle_output)
            print("\n appending shuffle output")

        final_training_output = pd.DataFrame(final_training_output, columns=['val_acc', 'val_auROC', 'plant', 'case', 'training size'])
        print("\n dataframe created \n")
        final_training_output.to_csv(f"results/{m_reads.split('_')[0]}_result.csv", index=False)