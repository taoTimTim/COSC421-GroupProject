# import mne
# import numpy as np
# import matplotlib.pyplot as plt
# from mne.connectivity import spectral_connectivity

# #from mne_connectivity import spectral_connectivity_epochs as spectral_connectivity



# # Load the BDF file
# raw = mne.io.read_raw_bdf("sub-001_task-med1breath_eeg.bdf", preload=True)

# raw.drop_channels(['EXG1', 'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8', 'GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp', 'Status'], on_missing='ignore')

# raw.filter(l_freq=0.5, h_freq=40.)


# # Extract data (channels x time)
# data = raw.get_data()
# ch_names = raw.info['ch_names']


# # Compute similarity (correlation) between electrodes
# corr_matrix = np.corrcoef(data)

# # Visualize correlation matrix
# plt.imshow(corr_matrix, cmap='viridis', interpolation='nearest')
# plt.colorbar(label='Correlation')
# plt.title('Electrode similarity (Pearson correlation)')
# plt.xticks(range(len(ch_names)), ch_names, rotation=90)
# plt.yticks(range(len(ch_names)), ch_names)
# plt.show()

# #print(ch_names)


# bands = {
#     'alpha': (8, 13),
#     'beta': (13, 30),
# }


# con_obj = spectral_connectivity(
# [data], method='wpli2_debiased', mode='fourier',
# sfreq=raw.info['sfreq'], fmin=8, fmax=13,
# faverage=True, verbose=False)

# # Extract dense matrix
# con = con_obj.get_data(output='dense')
# freqs = con_obj.freqs  # available if you need it

# print("Shape of wPLI matrix:", con.shape)

# # Plot the matrix
# plt.figure(figsize=(8, 6))
# plt.imshow(con, cmap='magma', vmin=0, vmax=1)
# plt.colorbar(label='wPLI')
# plt.title('Weighted Phase Lag Index (Alpha 8–13 Hz)')
# plt.xticks(range(len(ch_names)), ch_names, rotation=90)
# plt.yticks(range(len(ch_names)), ch_names)
# plt.tight_layout()
# plt.show()

import mne
import numpy as np
import matplotlib.pyplot as plt


from mne_connectivity import spectral_connectivity_epochs as spectral_connectivity

# ---- Load & clean ----
raw = mne.io.read_raw_bdf("data/sub-001_task-med1breath_eeg.bdf", preload=True)
raw.drop_channels(
    ['EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8',
     'GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp','Status'],
    on_missing='ignore'
)
raw.filter(l_freq=0.5, h_freq=40., verbose=False)
raw.resample(250, verbose=False)  # optional but helps speed

# ---- Make fixed-length epochs (best practice for spectral estimates) ----
# duration=2.0s with 50% overlap is a good default; tweak if you like
epochs = mne.make_fixed_length_epochs(raw, duration=2.0, overlap=1.0,preload=True, verbose=False)

# If you *really* want a single epoch spanning the whole recording, you can:
#epochs = mne.EpochsArray(raw.get_data()[None, ...], raw.info)  # 1 epoch
# But: fixed-length epochs are more stable for wPLI estimation.

# ---- Compute wPLI in alpha band ----
bands = {'alpha': (8, 13)}
fmin, fmax = bands['alpha']

# spectral_connectivity_epochs returns a SpectralConnectivity object
con_obj = spectral_connectivity(
    epochs,
    method='wpli2_debiased',
    mode='fourier',
    fmin=fmin, fmax=fmax,
    faverage=True,         # average across freqs in the band => 1 matrix
    mt_adaptive=False,     # multi-taper off is fine here
    n_jobs=1,              # or >1 if you want parallel
    verbose=False
)

# ---- Extract dense (channels x channels) matrix ----
con = con_obj.get_data(output='dense')  # 2D array
con = np.squeeze(con)                   # shape: (64, 64)
ch_names = epochs.info['ch_names']

print("wPLI matrix shape:", con.shape)   # (n_channels, n_channels)

# ---- Plot ----
plt.figure(figsize=(8, 6))
plt.imshow(con, cmap='magma', vmin=0, vmax=1)
plt.colorbar(label='wPLI')
plt.title('wPLI — Alpha (8–13 Hz)')
plt.xticks(range(len(ch_names)), ch_names, rotation=90)
plt.yticks(range(len(ch_names)), ch_names)
plt.tight_layout()
plt.show()

#---------------------------------------------------------------

import pandas as pd

# 'con' = your wPLI matrix (channels × channels)
# 'ch_names' = list of electrode names
channel_names = raw.info['ch_names']

# Create a DataFrame with labels
df_wpli = pd.DataFrame(con, index=channel_names, columns=channel_names)

# Option 1: Print full matrix (can be huge)
print(df_wpli)

# Option 2: Round values for readability
print(df_wpli.round(3))

# Option 3: Show only top-left part (if matrix is large)
print(df_wpli.round(3).iloc[:10, :10])
