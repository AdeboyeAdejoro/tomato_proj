import pandas as pd
import numpy as np
import shap
import tensorflow as tf
tf.compat.v1.disable_v2_behavior()  # Ensure compatibility with shap
from tensorflow.keras.models import load_model
from deeplift.visualization.viz_sequence import plot_weights_given_ax
from deeplift.visualization import viz_sequence
import matplotlib.pyplot as plt
import os
from utils import FastaSequenceLoader  # Adjust this import as per your utils module

# Ensure the results directory exists
results_path = 'results'
if not os.path.exists(results_path):
    os.mkdir(results_path)

# Define paths and variables
pennellii_genome = 'genomes/Spenn.fasta'
pennellii_gene_model = 'gene_models/spenn_v2.0_gene_models_annot.gff'
counts_file = "counts/pen_introgression_counts.csv"
pennellii_chromosomes = ['Spenn-ch01', 'Spenn-ch02', 'Spenn-ch03', 'Spenn-ch04', 'Spenn-ch05',
                  'Spenn-ch06', 'Spenn-ch07', 'Spenn-ch08', 'Spenn-ch09', 'Spenn-ch10', 'Spenn-ch11', 'Spenn-ch12']

# Create FastaSequenceLoader instance
fasta_loader = FastaSequenceLoader(pennellii_genome, pennellii_gene_model, val_chromosome=pennellii_chromosomes[0])

# Extract and encode sequences
encoded_train_seqs, encoded_val_seqs, train_ids, val_ids = fasta_loader.extract_seq()

# Combine training and validation sequences for prediction
#all_seqs = np.concatenate([encoded_train_seqs, encoded_val_seqs], axis=0)
#all_ids = train_ids + val_ids

val_sequences = np.stack(encoded_val_seqs)
print("Shape of val_sequences after stacking:", val_sequences.shape)




# Load the model
model_path = 'pen_model_Spenn-ch12_promoter_terminatorfruit_skin.h5'
model = load_model(model_path)

# Make predictions on the sequences
predictions = model.predict(val_sequences)

# Convert predictions to binary labels (assuming sigmoid activation in the output layer)
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
valid_indices = [i for i, gene_id in enumerate(val_ids) if gene_id in counts.index and counts.loc[gene_id, 'true_target'] != 2]
filtered_predictions = predicted_labels[valid_indices]
filtered_ids = [val_ids[i] for i in valid_indices]
filtered_true_targets = [counts.loc[gene_id, 'true_target'] for gene_id in filtered_ids]

# Calculate accuracy
accuracy = np.mean(np.array(filtered_predictions).flatten() == np.array(filtered_true_targets).flatten())
print(f'Accuracy: {accuracy}')


filters, biases = model.layers[0].get_weights()  # Adjust layer index as needed
print(f"Filter shape: {filters.shape}")

# Visualize filters
fig, ax = plt.subplots(2, 8, figsize=(30, 4))
ax_coords = [(i, j) for i in range(2) for j in range(8)]
for idx, ax_coord in enumerate(ax_coords):
    if idx < filters.shape[0]:
        plot_weights_given_ax(ax=ax[ax_coord], array=filters[idx],
                              height_padding_factor=0.2,
                              length_padding=1.0,
                              subticks_frequency=1,
                              highlight={})
    else:
        ax[ax_coord].axis('off')

plt.savefig(os.path.join(results_path, 'filter_visualization.png'))
plt.close(fig)

# Define a dinucleotide shuffling function
def dinuc_shuffle(onehot_seq, rng):
    indices = np.arange(onehot_seq.shape[0] // 2)
    rng.shuffle(indices)
    shuffled_seq = onehot_seq.copy()
    for i, j in zip(indices, indices[::-1]):
        shuffled_seq[2*i:2*i+2] = onehot_seq[2*j:2*j+2]
    return shuffled_seq

# Dinucleotide shuffling for multiple examples
def dinuc_shuffle_several_times(list_containing_input_modes_for_an_example, seed=1234):
    assert len(list_containing_input_modes_for_an_example) == 1
    onehot_seq = list_containing_input_modes_for_an_example[0]
    rng = np.random.RandomState(seed)
    to_return = np.array([dinuc_shuffle(onehot_seq, rng=rng) for _ in range(20)])
    return [to_return]

# Combine multiplier and difference from reference
def combine_mult_and_diffref(mult, orig_inp, bg_data):
    to_return = []
    for l in range(len(mult)):
        projected_hypothetical_contribs = np.zeros_like(bg_data[l]).astype("float")
        assert len(orig_inp[l].shape) == 2
        for i in range(orig_inp[l].shape[-1]):
            hypothetical_input = np.zeros_like(orig_inp[l]).astype("float")
            hypothetical_input[:, i] = 1.0
            hypothetical_difference_from_reference = (hypothetical_input[None, :, :] - bg_data[l])
            hypothetical_contribs = hypothetical_difference_from_reference * mult[l]
            projected_hypothetical_contribs[:, :, i] = np.sum(hypothetical_contribs, axis=-1)
        to_return.append(np.mean(projected_hypothetical_contribs, axis=0))
    return to_return
    
# Prepare sequences for explanation (positive and negative)
#pos_seqs = np.where(filtered_predictions == 1)[0]
#seqs_to_explain_1 = all_seqs[pos_seqs]  # All positive sequences
#print("seqs_to_explain_1:", seqs_to_explain_1.shape)

#neg_seqs = np.where(filtered_predictions == 0)[0]
#seqs_to_explain_2 = all_seqs[neg_seqs]  # All negative sequences
#print("seqs_to_explain_2:", seqs_to_explain_2.shape)



# Extract indices where predictions are correct
correct_pos_indices = [i for i in range(len(filtered_predictions)) if filtered_predictions[i] == 1 and filtered_true_targets[i] == 1]
correct_neg_indices = [i for i in range(len(filtered_predictions)) if filtered_predictions[i] == 0 and filtered_true_targets[i] == 0]

# Prepare sequences for explanation (positive and negative)
seqs_to_explain_1 = val_sequences[correct_pos_indices]
seqs_to_explain_2 = val_sequences[correct_neg_indices]

print("seqs_to_explain_1:", seqs_to_explain_1.shape)
print("seqs_to_explain_2:", seqs_to_explain_2.shape)



# Initialize the DeepExplainer
dinuc_shuff_explainer = shap.DeepExplainer(
    (model.input, model.output[:, 0]), 
    data=dinuc_shuffle_several_times, 
    combine_mult_and_diffref=combine_mult_and_diffref
)

# Compute SHAP values for positive sequences
shap_explanations_1 = dinuc_shuff_explainer.shap_values(seqs_to_explain_1)
print("shap_explanations_1 shape:", shap_explanations_1.shape)

# Multiply SHAP values by the original sequences and take mean
shap_explanations_1 = shap_explanations_1 * seqs_to_explain_1
shap_explanations_1 = shap_explanations_1.mean(axis=0)
print("shap_explanations_1 after mean operation:", shap_explanations_1.shape)

# Save the SHAP values visualization for positive sequences
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 2), sharey='row', sharex='col')
ax[0].plot(range(3020), shap_explanations_1[:, 0], c='cornflowerblue')
ax[0].plot(range(3020), shap_explanations_1[:, 1], c='green')
ax[0].plot(range(3020), shap_explanations_1[:, 2], c='crimson')
ax[0].plot(range(3020), shap_explanations_1[:, 3], c='orange')


# Compute SHAP values for negative sequences
shap_explanations_2 = dinuc_shuff_explainer.shap_values(seqs_to_explain_2)
print("shap_explanations_2 shape:", shap_explanations_2.shape)

# Multiply SHAP values by the original sequences and take mean
shap_explanations_2 = shap_explanations_2 * seqs_to_explain_2
shap_explanations_2 = shap_explanations_2.mean(axis=0)
print("shap_explanations_2 after mean operation:", shap_explanations_2.shape)

# Save the SHAP values visualization for negative sequences
ax[1].plot(range(3020), shap_explanations_2[:, 0], c='cornflowerblue')
ax[1].plot(range(3020), shap_explanations_2[:, 1], c='green')
ax[1].plot(range(3020), shap_explanations_2[:, 2], c='crimson')
ax[1].plot(range(3020), shap_explanations_2[:, 3], c='orange')
fig.tight_layout()
plt.savefig(os.path.join(results_path, 'saliency_map.png'))
plt.close(fig)
