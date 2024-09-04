import pandas as pd
import numpy as np
import shap
import tensorflow as tf
tf.compat.v1.disable_v2_behavior()  # Ensure compatibility with shap
from tensorflow.keras.models import load_model
from deeplift.visualization.viz_sequence import plot_weights_given_ax
from deeplift.visualization import viz_sequence
import h5py
import matplotlib.pyplot as plt
import os
from utils import FastaSequenceLoader  
import modisco
from tensorflow.keras import backend
from importlib import reload

os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
os.environ['TF_GPU_ALLOCATOR'] = 'cuda_malloc_async'


# Ensure the results directory exists
results_path = 'modisco_results'
if not os.path.exists(results_path):
    os.mkdir(results_path)

# Define paths and variables
pennellii_genome = 'genomes/S_lycopersicum_chromosomes.4.00.fa'
pennellii_gene_model = 'gene_models/ITAG4.0_gene_models.gff'
counts_file = "counts/pure_lyc_counts.csv"
pennellii_chromosomes = ['SL4.0ch01', 'SL4.0ch02', 'SL4.0ch03', 'SL4.0ch04', 'SL4.0ch05',
                  'SL4.0ch06', 'SL4.0ch07', 'SL4.0ch08', 'SL4.0ch09', 'SL4.0ch10', 'SL4.0ch11', 'SL4.0ch12']

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
model_path = 'pure_model_SL4.0ch12_promoter_terminatorfruit_skin.h5'
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
shap_explanations_1 = dinuc_shuff_explainer.shap_values(seqs_to_explain_1) #hypo
print("shap_explanations_1 shape:", shap_explanations_1.shape)

# Multiply SHAP values by the original sequences and take mean
actual_shap_explanations_1 = shap_explanations_1 * seqs_to_explain_1   #contrib
#actual_shap_explanations_1 = shap_explanations_1.mean(axis=0)
print("actual_shap_explanations_1 without mean operation:", actual_shap_explanations_1.shape)


# Compute SHAP values for negative sequences
shap_explanations_2 = dinuc_shuff_explainer.shap_values(seqs_to_explain_2)
print("shap_explanations_2 shape:", shap_explanations_2.shape)

# Multiply SHAP values by the original sequences and take mean
actual_shap_explanations_2 = shap_explanations_2 * seqs_to_explain_2
actual_shap_explanations_2 = shap_explanations_2.mean(axis=0)
print("actual_shap_explanations_2 after mean operation:", actual_shap_explanations_2.shape)

def run_modisco(species):
    save_file = f"modisco_results/{species}_modisco.hdf5"

    contribution_scores = actual_shap_explanations_1
    hypothetical_scores = shap_explanations_1
    one_hots = seqs_to_explain_1

    print('contributions', contribution_scores.shape)
    print('hypothetical scores', hypothetical_scores.shape)
    print('correct predictions', one_hots.shape)


    reload(modisco.util)
    reload(modisco.pattern_filterer)
    reload(modisco.aggregator)
    reload(modisco.core)
    reload(modisco.seqlet_embedding.advanced_gapped_kmer)
    reload(modisco.affinitymat.transformers)
    reload(modisco.affinitymat.core)
    reload(modisco.affinitymat)
    reload(modisco.cluster.core)
    reload(modisco.cluster)
    reload(modisco.tfmodisco_workflow.seqlets_to_patterns)
    reload(modisco.tfmodisco_workflow)
    reload(modisco)

    #reshaped_contribution_scores = np.tile(contribution_scores, (hypothetical_scores.shape[0], 1, 1))

    null_per_pos_scores = modisco.coordproducers.LaplaceNullDist(num_to_samp=5000)
    tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
        sliding_window_size = 21,
        flank_size = 5,
        target_seqlet_fdr=0.01,
        max_seqlets_per_metacluster=50000,
        seqlets_to_patterns_factory=modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
            trim_to_window_size=10,
            initial_flank_to_add=2,
            final_flank_to_add=0,
            final_min_cluster_size=60,
            n_cores=50)

    )(
        task_names=['task0'],
        contrib_scores={'task0': contribution_scores},
        hypothetical_contribs={'task0': hypothetical_scores},
        one_hot=one_hots,
        null_per_pos_scores=null_per_pos_scores
     )

    reload(modisco.util)
    grp = h5py.File(save_file, "w")
    tfmodisco_results.save_hdf5(grp)
    grp.close()



species = ["lycopersicum"]
genomes = ['genomes/S_lycopersicum_chromosomes.4.00.fa']
gene_models = ['gene_models/ITAG4.0_gene_models.gff'] 
mapped_read_counts = ["counts/pure_lyc_counts.csv"]


for plant, fasta_file, gtf_file, counts in zip(species, genomes, gene_models, mapped_read_counts):
    if not os.path.exists(f'modisco_results/{plant}_modisco.hdf5'):
        run_modisco(plant)
