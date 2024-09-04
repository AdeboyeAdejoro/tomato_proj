import os
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from modisco.visualization.viz_sequence import plot_weights_given_ax

if not os.path.exists('figures/motifs'):
    os.mkdir('figures/motifs')

def get_predictive_pwms(mod_file: str, species: str):
    print(species)
    cwm = []  # Ensure this is a list
    motif_id = []
    strand = []
    metacluster_id = []
    n_motif_seqlets = []

    f = h5py.File(mod_file, 'r')
    for metacluster_idx, metacluster_key in enumerate(f["metacluster_idx_to_submetacluster_results"]):
        print(metacluster_idx, metacluster_key)
        metacluster = f["metacluster_idx_to_submetacluster_results"][metacluster_key]
        patterns = metacluster["seqlets_to_patterns_result"]["patterns"]

        for pattern_idx, pattern_name in enumerate(patterns['all_pattern_names']):
            pattern = patterns[pattern_name.decode()]
            pattern_seqlets = pattern['seqlets_and_alnmts']['seqlets']

            # Forward strand
            nfcwm = np.absolute(pattern["task0_contrib_scores"]["fwd"][:])
            if nfcwm.size > 0:  # Ensure non-empty arrays
                nfcwm = len(pattern_seqlets) * (nfcwm / np.max(nfcwm.flat))
                cwm.append(nfcwm)
                motif_id.append(pattern_name.decode())
                metacluster_id.append(metacluster_key)
                n_motif_seqlets.append(len(pattern_seqlets))
                strand.append('fwd')
                save_fwd = f'figures/motifs/{species}/{metacluster_key}_{pattern_name.decode()}_fwd.png'
                plot_and_save_weights(pattern["task0_contrib_scores"]["fwd"][:], save_fwd)

            # Reverse strand
            if "rev" in pattern["task0_contrib_scores"]:
                nrcwm = np.absolute(pattern["task0_contrib_scores"]["rev"][:])
                if nrcwm.size > 0:  # Ensure non-empty arrays
                    nrcwm = len(pattern_seqlets) * (nrcwm / np.max(nrcwm.flat))
                    cwm.append(nrcwm)
                    motif_id.append(pattern_name.decode())
                    metacluster_id.append(metacluster_key)
                    n_motif_seqlets.append(len(pattern_seqlets))
                    strand.append('rev')
                    save_rev = f'figures/motifs/{species}/{metacluster_key}_{pattern_name.decode()}_rev.png'
                    plot_and_save_weights(pattern["task0_contrib_scores"]["rev"][:], save_rev)

    cwm = np.array(cwm)  # Convert to NumPy array only once all data is appended
    motif_save = f'figures/motifs/{species}/motifs.h5'
    with h5py.File(motif_save, 'w') as h:
        h.create_dataset('CWMs', data=cwm)

    meta_info = pd.DataFrame({
        'motifID': motif_id,
        'strand': strand,
        'metacluster': metacluster_id,
        'n_seqlets': n_motif_seqlets
    })
    meta_info.to_csv(f'figures/motifs/{species}/meta_info.csv', sep="\t", index=False)

def plot_and_save_weights(weights: np.ndarray, save_path: str):
    # Ensure the weights are in the correct shape and type
    if weights.ndim != 2 or not np.isfinite(weights).all():
        raise ValueError(f"Invalid weights shape or non-finite values: shape {weights.shape}, finite {np.isfinite(weights).all()}")

    # Plot weights and save the figure
    fig = plt.figure(figsize=(10, 2))
    ax = fig.gca()
    plot_weights_given_ax(ax=ax, array=weights)
    plt.savefig(save_path)
    plt.close(fig)

modisco_feats = ['lycopersicum_modisco.hdf5']
species = ["Lyc"]

for feats, sp in zip(modisco_feats, species):
    if not os.path.exists(f'figures/motifs/{sp}'):
        os.mkdir(f'figures/motifs/{sp}')
    feats_path = f'modisco_results/{feats}'
    get_predictive_pwms(feats_path, sp)
