import mne
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from mne_connectivity import spectral_connectivity_epochs as spectral_connectivity


# Configuration
data_path = "./data"
results_path = "./results"
os.makedirs(results_path, exist_ok=True)

# EEG processing parameters
l_freq = 0.5
h_freq = 40.
resample_sfreq = 250
epoch_duration = 2.0
epoch_overlap = 1.0
bands = {'alpha': (8, 13), 'beta': (13, 30)}
fmin, fmax = bands['alpha']

# Helper function to clean raw EEG
def preprocess_raw(raw):
    # Drop non-EEG channels
    raw.drop_channels(
        ['EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8',
         'GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp','Status'],
        on_missing='ignore'
    )
    raw.filter(l_freq=l_freq, h_freq=h_freq, verbose=False)
    raw.resample(resample_sfreq, verbose=False)
    return raw

# Process each file 
for fname in os.listdir(data_path):
    if not fname.endswith(".bdf"):
        continue

    print(f"\nProcessing {fname}...")

    # Load raw EEG
    raw = mne.io.read_raw_bdf(os.path.join(data_path, fname), preload=True)
    raw = preprocess_raw(raw)

    # Make epochs
    epochs = mne.make_fixed_length_epochs(
        raw, duration=epoch_duration, overlap=epoch_overlap, preload=True, verbose=False
    )

    # Determine task from filename
    if "task-med1breath" in fname:
        task_name = "med1breath"
    elif "task-med2" in fname:
        task_name = "med2"
    elif "task-think" in fname:
        task_name = "thinking"
    else:
        task_name = "unknown"

    # Loop over frequency bands
    for band_name, (fmin, fmax) in bands.items():

        print(f"  Computing {band_name} band ({fmin}-{fmax} Hz)...")

        # Compute wPLI
        con_obj = spectral_connectivity(
            epochs,
            method='wpli2_debiased',
            mode='fourier',
            fmin=fmin, fmax=fmax,
            faverage=True,
            mt_adaptive=False,
            n_jobs=1,
            verbose=False
        )

        # Extract dense symmetric matrix
        con = np.squeeze(con_obj.get_data(output='dense'))
        con = (con + con.T) / 2

        # Create save folder: e.g., results/alpha_med1breath/
        save_folder = os.path.join(results_path, f"{band_name}_{task_name}")
        os.makedirs(save_folder, exist_ok=True)

        # Save CSV
        channel_names = epochs.ch_names
        df_wpli = pd.DataFrame(con, index=channel_names, columns=channel_names)
        csv_name = fname.replace(".bdf", f"_{band_name}_wpli.csv")
        df_wpli.to_csv(os.path.join(save_folder, csv_name), index=False)

        print(f"    Saved {csv_name} to {save_folder}")


#Load & clean

# subject_file = "sub-001_task-med1breath_eeg.bdf"

# raw = mne.io.read_raw_bdf(subject_file, preload=True)
# # Drop non-EEG channels
# raw.drop_channels(
#     ['EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8',
#      'GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp','Status'],
#     on_missing='ignore'
# )
# raw.filter(l_freq=0.5, h_freq=40., verbose=False)
# raw.resample(250, verbose=False) 

# # Make fixed-length epochs (duration=2.0s with 50% overlap)
# epochs = mne.make_fixed_length_epochs(raw, duration=2.0, overlap=1.0,preload=True, verbose=False)


# # ---- Compute wPLI in alpha band ----
# bands = {'alpha': (8, 13)}
# fmin, fmax = bands['alpha']



# # spectral_connectivity_epochs returns a SpectralConnectivity object
# con_obj = spectral_connectivity(
#     epochs,
#     method='wpli2_debiased',
#     mode='fourier',
#     fmin=fmin, fmax=fmax,
#     faverage=True,        
#     mt_adaptive=False,    
#     n_jobs=1,            
#     verbose=False
# )

# #Extract dense (channels x channels) matrix
# con = con_obj.get_data(output='dense')  # 2D array
# con = np.squeeze(con)                   # shape: (64, 64)
# con = (con + con.T) / 2
# print("Symmetric?", np.allclose(con, con.T, atol=1e-6))

# channel_names = epochs.ch_names




# # ---- Plot ----
# # plt.figure(figsize=(8, 6))
# # plt.imshow(con, cmap='magma', vmin=0, vmax=1)
# # plt.colorbar(label='wPLI')
# # plt.title('wPLI — Alpha (8–13 Hz)')
# # plt.xticks(range(len(ch_names)), ch_names, rotation=90)
# # plt.yticks(range(len(ch_names)), ch_names)
# # plt.tight_layout()
# # plt.show()

# #---------------------------------------------------------------



# output_path = "./results"

# os.makedirs(output_path, exist_ok=True)

# # 'con' = your wPLI matrix (channels × channels)
# # 'ch_names' = list of electrode names


# # Create a DataFrame with labels
# df_wpli = pd.DataFrame(con, index=channel_names, columns=channel_names)

# # Save to CSV
# csv_name = subject_file.replace(".bdf", "_wpli.csv")
# df_wpli.to_csv(os.path.join(output_path, csv_name))

# #Print full matrix
# #print(df_wpli)
