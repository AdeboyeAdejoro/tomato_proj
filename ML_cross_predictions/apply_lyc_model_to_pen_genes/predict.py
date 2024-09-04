import pandas as pd
import numpy as np
from utils import FastaSequenceLoader
from tensorflow.keras.models import load_model

# Define paths to the Pennellii genome and gene model files
pennellii_genome = 'genomes/Spenn.fasta'
pennellii_gene_model = 'gene_models/spenn_v2.0_gene_models_annot.gff'
counts_file = "counts/pen_introgression_counts.csv"

# Define the chromosomes of interest (example chromosomes)
pennellii_chromosomes = ['Spenn-ch01', 'Spenn-ch02', 'Spenn-ch03', 'Spenn-ch04', 'Spenn-ch05',
                  'Spenn-ch06', 'Spenn-ch07', 'Spenn-ch08', 'Spenn-ch09', 'Spenn-ch10', 'Spenn-ch11', 'Spenn-ch12']

# Create FastaSequenceLoader instance
fasta_loader = FastaSequenceLoader(pennellii_genome, pennellii_gene_model, val_chromosome=pennellii_chromosomes[0])  # Replace with appropriate chromosome

# Extract and encode sequences
encoded_train_seqs, encoded_val_seqs, train_ids, val_ids = fasta_loader.extract_seq()

# Combine training and validation sequences for prediction
all_seqs = np.array(encoded_train_seqs + encoded_val_seqs)
all_ids = train_ids + val_ids

# Path to the best performing model
best_model_path = 'pure_model_SL4.0ch12_promoter_terminatorfruit_skin.h5'

# Load the model
best_model = load_model(best_model_path)

# Make predictions on the Pennellii sequences
predictions = best_model.predict(all_seqs)

# Convert predictions to binary labels (if using sigmoid activation in the output layer)
predicted_labels = (predictions > 0.5).astype(int)


# Load the counts file and generate true labels
counts = pd.read_csv(counts_file, index_col=0)
true_targets = []

for log_count in counts['median'].values:
    if log_count <= np.percentile(counts['median'], 25):
        true_targets.append(0)
    elif log_count >= np.percentile(counts['median'], 75):
        true_targets.append(1)
    else:
        true_targets.append(2)

counts['true_target'] = true_targets




# Filter out the predictions and true targets to only include genes with defined labels
valid_indices = [i for i, gene_id in enumerate(all_ids) if gene_id in counts.index and counts.loc[gene_id, 'true_target'] != 2]
filtered_predictions = predicted_labels[valid_indices]
filtered_ids = [all_ids[i] for i in valid_indices]
filtered_true_targets = [counts.loc[gene_id, 'true_target'] for gene_id in filtered_ids]

# Calculate accuracy
accuracy = np.mean(np.array(filtered_predictions).flatten() == np.array(filtered_true_targets).flatten())
print(f'Accuracy: {accuracy}')

# Create a DataFrame to save the predictions and true labels
predictions_df = pd.DataFrame({
    'gene_id': filtered_ids,
    'prediction': filtered_predictions.flatten(),
    'probability': predictions[valid_indices].flatten(),
    'true_label': filtered_true_targets
})

# Save the predictions to a CSV file
predictions_df.to_csv('results/pennellii_predictions_with_labels.csv', index=False)

print(predictions_df.head())