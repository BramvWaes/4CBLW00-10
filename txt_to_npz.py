import os
import numpy as np

def parse_ir_file(filepath):
    """
    Parses a single .txt IR spectrum file.

    Args:
        filepath (str): The path to the .txt file.

    Returns:
        tuple: A tuple containing:
            - labels (dict): A dictionary of chemical group labels (e.g., {'Alkene': 0, ...}).
            - datapoints (numpy.ndarray): A 2D NumPy array of [X, Y] coordinates.
    """
    labels = {}
    datapoints = []
    in_compounds_section = False
    in_datapoints_section = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line == '#START COMPOUNDS':
                in_compounds_section = True
                in_datapoints_section = False
                continue
            elif line == '#END COMPOUNDS':
                in_compounds_section = False
                continue
            elif line == '#START DATAPOINTS':
                in_datapoints_section = True
                in_compounds_section = False
                continue
            elif line == '#END DATAPOINTS':
                in_datapoints_section = False
                continue

            if in_compounds_section:
                try:
                    group, value = line.split(':')
                    labels[group.strip()] = int(value.strip())
                except ValueError:
                    print(f"Warning: Could not parse label line '{line}' in {filepath}")
            elif in_datapoints_section:
                try:
                    x, y = map(float, line.split())
                    datapoints.append([x, y])
                except ValueError:
                    print(f"Warning: Could not parse datapoint line '{line}' in {filepath}")

    return labels, np.array(datapoints)

def convert_folder_to_npz(input_folder, output_filepath):
    """
    Converts all .txt IR spectrum files in a folder to a single .npz file.

    Args:
        input_folder (str): Path to the folder containing the .txt sample files.
        output_filepath (str): Path to save the output .npz file.
    """
    all_data = []
    all_labels = []
    all_filenames = [] # To keep track of original filenames

    # Get a consistent order of chemical groups from the first file
    first_file_path = None
    for filename in os.listdir(input_folder):
        if filename.endswith('.txt'):
            first_file_path = os.path.join(input_folder, filename)
            break

    if not first_file_path:
        print(f"No .txt files found in {input_folder}. Exiting.")
        return

    # Parse the first file to get the ordered list of all possible labels
    sample_labels, _ = parse_ir_file(first_file_path)
    # Get sorted keys to ensure consistent label order across all samples
    label_keys = sorted(sample_labels.keys())
    print(f"Identified chemical groups (labels): {label_keys}")

    print(f"Processing files in: {input_folder}")
    for filename in os.listdir(input_folder):
        if filename.endswith('.txt'):
            filepath = os.path.join(input_folder, filename)
            try:
                labels, datapoints = parse_ir_file(filepath)

                # Convert labels dictionary to a consistent numpy array based on label_keys
                current_labels_array = np.array([labels.get(key, 0) for key in label_keys], dtype=np.int32)

                all_data.append(datapoints)
                all_labels.append(current_labels_array)
                all_filenames.append(filename)
                print(f"Processed {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    if not all_data:
        print("No data was successfully processed. NPZ file not created.")
        return

    # Pad data if all_data arrays have different lengths
    # find the max length and pad with zeros.

    max_len = max(data.shape[0] for data in all_data)
    padded_data = np.zeros((len(all_data), max_len, 2), dtype=np.float32) # 2 for X, Y coordinates

    for i, data in enumerate(all_data):
        padded_data[i, :data.shape[0], :] = data

    # Convert lists to numpy arrays for saving
    final_labels = np.array(all_labels, dtype=np.int32)

    # Save to .npz file
    np.savez_compressed(output_filepath,
                         data=padded_data,
                         labels=final_labels,
                         label_names=np.array(label_keys), # Save the names of the labels for reference
                         filenames=np.array(all_filenames))

    print(f"\nSuccessfully converted {len(all_data)} samples to {output_filepath}")
    print(f"Data shape: {padded_data.shape}")
    print(f"Labels shape: {final_labels.shape}")

if __name__ == "__main__":
    folder = "TRANSMISSION" #foldet containing the .txt files
    output_npz_file = "ir_spectra_dataset.npz" #Output file name

    # Run the conversion
    convert_folder_to_npz(folder, output_npz_file)

    #Verifying the file manually 
    print("\n--- Verifying the generated NPZ file ---")
    try:
        loaded_data = np.load(output_npz_file)
        print(f"Keys in NPZ file: {list(loaded_data.keys())}")
        print(f"Loaded data shape: {loaded_data['data'].shape}")
        print(f"Loaded labels shape: {loaded_data['labels'].shape}")
        print(f"Loaded label names: {loaded_data['label_names']}")
        print(f"First sample data (first 5 rows):\n{loaded_data['data'][0, :5, :]}")
        print(f"First sample labels:\n{loaded_data['labels'][0]}")
        print(f"Second sample data (first 5 rows):\n{loaded_data['data'][1, :5, :]}")
        print(f"Second sample labels:\n{loaded_data['labels'][1]}")
        print(f"Third sample data (first 5 rows):\n{loaded_data['data'][2, :5, :]}")
        print(f"Third sample labels:\n{loaded_data['labels'][2]}")

    except Exception as e:
        print(f"Error loading NPZ file for verification: {e}")

