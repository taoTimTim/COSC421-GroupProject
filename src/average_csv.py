import os
import pandas as pd
import numpy as np

results_path = "./results"
averages_path = os.path.join(results_path, "averages")
os.makedirs(averages_path, exist_ok=True)

# Loop through all subfolders in results/
for subfolder in os.listdir(results_path):
    folder_path = os.path.join(results_path, subfolder)

    # Skip files and skip the main averages folder itself
    if not os.path.isdir(folder_path) or subfolder == "averages":
        continue

    print(f"\nProcessing folder: {subfolder}")

    # Collect CSV files
    csv_files = [f for f in os.listdir(folder_path) if f.endswith(".csv")]
    if not csv_files:
        print("  No CSV files found, skipping...")
        continue

    matrices = []
    channel_index = None
    channel_cols = None

    for csv_file in csv_files:
        df = pd.read_csv(os.path.join(folder_path, csv_file), index_col=0)
        matrices.append(df.values)
        channel_index = df.index
        channel_cols = df.columns

    matrices = np.array(matrices)
    avg_matrix = np.mean(matrices, axis=0)

    # Create subfolder in averages/
    avg_subfolder = os.path.join(averages_path, subfolder)
    os.makedirs(avg_subfolder, exist_ok=True)

    #ADD THE NUMBER RANGE OF SUBJECTS IN THE FILE NAME
    sub_range = "sub1-2"  # MODIFY THIS BASED ON ACTUAL SUBJECTS INCLUDED

    # Save file
    avg_df = pd.DataFrame(avg_matrix, index=channel_index, columns=channel_cols)
    avg_csv_name = os.path.join(avg_subfolder, f"{subfolder}_wpli_average_{sub_range}.csv")
    avg_df.to_csv(avg_csv_name)

    print(f" Saved averaged CSV → {avg_csv_name}")
