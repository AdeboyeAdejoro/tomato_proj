import os
import pandas as pd
import numpy as np
from pyfaidx import Fasta
import pyranges as pr
from tensorflow.keras import Sequential, optimizers, backend, models
from tensorflow.keras.layers import Conv1D, Dense, MaxPool1D, Dropout, Flatten
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.utils import shuffle


def onehot(seq):
    code = {'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'unk':[0, 0, 0, 0]}
    encoded = np.zeros(shape=(len(seq), 4))
    for i, nt in enumerate(seq):
        if nt in ['A', 'C', 'G', 'T']:
            encoded[i, :] = code[nt]
        else:
            encoded[i, :] = code['unk']
    return encoded

class FastaSequenceLoader:
    def __init__(self, fasta, gff, val_chromosome, upstream=1000, downstream=500):
        self.fasta = Fasta(fasta, as_raw=True, sequence_always_upper=True, read_ahead=10000)
        gene_models = pr.read_gff3(gff, as_df=True)
        gene_models = gene_models[gene_models['Feature'] == 'gene']
        gene_models['ID'] = gene_models['ID'].str.replace('gene:', '')
        gene_models = gene_models[['Chromosome', 'Start', 'End', 'Strand', 'ID']]
        self.gff = gene_models
        self.val_chromosome = val_chromosome
        self.upstream = upstream
        self.downstream = downstream
        
    
    def extract_seq(self):
        print(self.gff.head())
        encoded_val_seqs, val_ids  = [], []
        encoded_train_seqs, train_ids = [], []
        for chrom, start, end, strand, gene_id in self.gff.values:
            if strand == "+":
                prom_start, prom_end = start - self.upstream, start + self.downstream
                term_start, term_end = end - self.downstream, end + self.upstream
                if prom_start > 0 and term_start > 0:
                    encoded_seq = np.concatenate([onehot(self.fasta[chrom][prom_start:prom_end]), np.zeros(shape=(20, 4)), onehot(self.fasta[chrom][term_start:term_end])])
                    if encoded_seq.shape[0] == 2 * (self.upstream + self.downstream) + 20:
                        if chrom == self.val_chromosome:
                            encoded_val_seqs.append(encoded_seq)
                            val_ids.append(gene_id)
                            
                        else:
                            encoded_train_seqs.append(encoded_seq)
                            train_ids.append(gene_id)
                            
                            
            else:
                prom_start, prom_end = end - self.downstream, end + self.upstream
                term_start, term_end = start - self.upstream, start + self.downstream
                if prom_start > 0 and term_start > 0:
                    encoded_seq = np.concatenate([onehot(self.fasta[chrom][prom_start:prom_end])[::-1, ::-1], np.zeros(shape=(20, 4)), onehot(self.fasta[chrom][term_start:term_end])[::-1, ::-1]])
                    if encoded_seq.shape[0] == 2 * (self.upstream + self.downstream) + 20:
                        if chrom == self.val_chromosome:
                            encoded_val_seqs.append(encoded_seq)
                            val_ids.append(gene_id)
                            
                        else:
                            encoded_train_seqs.append(encoded_seq)
                            train_ids.append(gene_id)
                            
                            
                            
        return encoded_train_seqs, encoded_val_seqs, train_ids, val_ids
                            
                            

    
    
    
class ConvNetwork:
    def __init__(self, encoded_train_seqs, encoded_val_seqs, train_ids, val_ids, val_chromosome, counts, organism, case, outer_flank=1000, inner_flank=500, tissue=""):
        self.train_seqs = encoded_train_seqs
        self.val_seqs = encoded_val_seqs
        self.train_ids = train_ids
        self.val_ids = val_ids
        self.val_chrom = val_chromosome
        self.counts = counts
        self.organism = organism
        self.case = case
        self.outer = outer_flank
        self.inner = inner_flank
        self.tissue = tissue
        
    
    
    def build_network(self, x_train, x_val, y_train, y_val):
        backend.clear_session()
        model = Sequential([
            #Conv Block 1
            Conv1D(64, kernel_size=8, activation='relu', padding='same', input_shape=(x_train.shape[1], x_train.shape[2])),
            Conv1D(64, kernel_size=8, activation='relu', padding='same'),
            MaxPool1D(8, padding='same'),
            Dropout(0.25),
            
            #Conv Block 2
            Conv1D(128, kernel_size=8, activation='relu', padding='same'),
            Conv1D(128, kernel_size=8, activation='relu', padding='same'),
            MaxPool1D(8, padding='same'),
            Dropout(0.25),
            
            #Conv Block 3
            
            Conv1D(64, kernel_size=8, activation='relu', padding='same'),
            Conv1D(64, kernel_size=8, activation='relu', padding='same'),
            MaxPool1D(8, padding='same'),
            Dropout(0.25),
            
            #Fully connected Block
            Flatten(),
            Dense(128, activation='relu'),
            Dropout(0.25),
            Dense(64, activation='relu'),
            Dense(1, activation='sigmoid')
            
        ])
        
        print(model.summary())
        model_save_name = f"saved_models/{self.organism}_model_{self.val_chrom}_{self.case}_{self.tissue}.h5"
        print(f"save path: {model_save_name}")
        model_chkpt = ModelCheckpoint(model_save_name, save_best_only=True, verbose=1)
        early_stop = EarlyStopping(patience=10)
        reduce_lr = ReduceLROnPlateau(patience=5, factor=0.1)
        model.compile(loss='binary_crossentropy', optimizer=optimizers.Adam(0.0001), metrics=['accuracy'])
        model.fit(x_train, y_train, batch_size=64, epochs=100, validation_data=(x_val, y_val), callbacks=[early_stop, model_chkpt, reduce_lr])
        
        saved_model = models.load_model(model_save_name)
        predictions = saved_model.predict(x_val)
        val_auroc = roc_auc_score(y_val, predictions)
        predictions = predictions > 0.5
        val_acc = accuracy_score(y_val, predictions)
        print('Best model performance-------------\n')
        print(f'Accuracy: {val_acc}, auROC: {val_auroc}\n')
        print('----------------------------------------')
        
        performance = [val_acc, val_auroc, self.organism, self.case, x_train.shape[0]]
        return performance
    
    def train_network(self):
        train_labels, train_seqs = [], []
        val_labels, val_seqs = [], []
        for train_id, train_seq in zip(self.train_ids, self.train_seqs):
            train_labels.append(self.counts.loc[train_id, 'true_target'])
            train_seqs.append(train_seq)
            
        for val_id, val_seq in zip(self.val_ids, self.val_seqs):
            val_labels.append(self.counts.loc[val_id, 'true_target'])
            val_seqs.append(val_seq)
            
        train_labels, val_labels = np.array(train_labels), np.array(val_labels)
        train_seqs, val_seqs, = np.array(train_seqs), np.array(val_seqs)
        
        #Random down sampling to balance data
        low_train, high_train = np.where(train_labels== 0)[0], np.where(train_labels == 1)[0]
        print(f"low_train: {len(low_train)}, high_train: {len(high_train)}")
        min_class = min([len(low_train), len(high_train)])
        selected_low_train = np.random.choice(low_train, min_class, replace=False)
        selected_high_train = np.random.choice(high_train, min_class, replace=False)
        x_train = np.concatenate([np.take(train_seqs, selected_low_train, axis=0), np.take(train_seqs, selected_high_train, axis=0)], axis=0)
        y_train = np.concatenate([np.take(train_labels, selected_low_train, axis=0), np.take(train_labels, selected_high_train, axis=0)], axis=0)
        x_train, y_train = shuffle(x_train, y_train, random_state=42)
        
        #selecting validation sequences with label 1 and 0
        low_val, high_val = np.where(val_labels ==0)[0], np.where(val_labels == 1)[0]
        print(f"low_val: {len(low_val)}, high_val: {len(high_val)}")
        x_val = np.concatenate([np.take(val_seqs, low_val, axis=0), np.take(val_seqs, high_val, axis=0)])
        y_val = np.concatenate([np.take(val_labels, low_val, axis=0), np.take(val_labels, high_val, axis=0)])
        print(x_train.shape, x_val.shape)
        print(f"validation size: {x_val.shape[0]}")
        
        #Masking first three nucleotides after gene start and before gene end
        x_train[:, self.outer:self.outer +3, :] = 0
        x_train[:, self.outer + (self.inner * 2) + 17 : self.outer + (self.inner * 2) + 20, :] = 0
        x_val[:, self.outer:self.outer + 3, :] = 0
        x_val[:, self.outer + (self.inner * 2) + 17: self.outer + (self.inner * 2) + 20, :] = 0
        
        output = self.build_network(x_train, x_val, y_train, y_val)
        return output
    
    
    
    
    

    
  
