from jcamp import jcamp_readfile, jcamp_writefile
from os import listdir
from pathlib import Path
import os
import numpy as np
import concurrent.futures as cf
import json
import logging
from rich.console import Console
from rich.progress import Progress

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[logging.StreamHandler()])
logger = logging.getLogger(__name__)
def convert_to_wavenumbers(value):
    """
    Convert a single value from micrometers to wavenumbers (1/cm).
    """
    return 1e4 / value

def write_to_txt(datapoints, groups, filename):
    output_text = ""
    output_text += "#START COMPOUNDS\n"
    for compound, present in groups.items():
        output_text += f"{compound}: {present}\n"
    output_text += "#END COMPOUNDS\n#START DATAPOINTS\n"
    for i in range(len(datapoints["x"])):
        x = datapoints["x"][i]
        y = datapoints["y"][i]
        output_text += f"{x} {y}\n"
    output_text += "#END DATAPOINTS\n"
    with open(filename, "w") as filehandle:
        filehandle.write(output_text)
    logger.debug(f"Data written to {filename}")


def transmittance_to_absorbance(y_vals):
    """
    Convert transmittance values to absorbance.

    Args:
        y_vals (list or np.array): Transmittance values (0–1 or 0–100).

    Returns:
        list: Absorbance values.
    """
    y = np.array(y_vals, dtype=np.float64)

    # Normalize to range [0–1] if in percent
    if np.max(y) > 1.0:
        y = y / 100.0

    # Avoid log(0) by setting a small floor value
    y = np.clip(y, 1e-6, 1.0)

    absorbance = -np.log10(y)
    return absorbance.tolist()

def get_datapoint_count(file_path):
    """
    Get the number of data points in a JCAMP file.

    Args:
        file_path (str): Path to the JCAMP file.

    Returns:
        int: Number of data points.
    """
    jcamp_data = jcamp_readfile(file_path)
    return jcamp_data["npoints"]

def normalize_file(file_path, groups):
    logger.debug(f"Normalizing file: {file_path}")
    jcamp_data = jcamp_readfile(file_path)
    sampling_procedure = str(jcamp_data.get("sampling procedure", "Unknown"))

    if jcamp_data['xunits'] == "MICROMETERS":
        # Convert x values to wavenumbers
        new_x = np.array([convert_to_wavenumbers(x) for x in jcamp_data["x"]])
        jcamp_data["x"] = list(reversed(new_x.tolist()))
        jcamp_data["xunits"] = "1/CM"

        # Update metadata
        if jcamp_data["xfactor"] != 1:
            print(f"Warning: xfactor is not 1 for {file}. It will be set to 1.")
        jcamp_data["xfactor"] = 1  # Assuming 1 for simplicity; adjust if needed
        jcamp_data["deltax"] = float(np.mean(np.diff(new_x)))
        jcamp_data["firstx"] = float(new_x[0])
        jcamp_data["lastx"] = float(new_x[-1])
        jcamp_data["minx"] = float(np.min(new_x))
        jcamp_data["maxx"] = float(np.max(new_x))
    if jcamp_data["yunits"] == "TRANSMITTANCE":
        # Convert y values to absorbance
        jcamp_data["y"] = transmittance_to_absorbance(jcamp_data["y"])
        jcamp_data["yunits"] = "ABSORBANCE"
        jcamp_data["ymin"] = float(np.min(jcamp_data["y"]))
        jcamp_data["ymax"] = float(np.max(jcamp_data["y"]))
        jcamp_data["miny"] = float(np.min(jcamp_data["y"]))
        jcamp_data["maxy"] = float(np.max(jcamp_data["y"]))
        jcamp_data["deltay"] = float(np.mean(np.diff(jcamp_data["y"])))
        jcamp_data["firsty"] = float(jcamp_data["y"][0])
        jcamp_data["lasty"] = float(jcamp_data["y"][-1])
        jcamp_data["yfactor"] = 1.0

    if jcamp_data["yunits"] == "ABSORBANCE":
        # Write back to a categorized directory
        file_name = Path(file_path).stem
        output_dir = os.path.join(folder_path, sampling_procedure)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{file_name}.txt")
        logger.debug(f"Writing normalized data to {output_path}")
        write_to_txt(jcamp_data, groups, output_path)
def normalize_folder(folder_path, groups):
    """
    Normalize all JCAMP files in the specified folder and update metadata fields.

    Args:
        folder_path (str): Path to the folder containing JDX files.
    """
    # List all files in the folder
    files = listdir(folder_path)
    
    logger.debug(f"Starting normalization for {len(files)} files in {folder_path}")
    with Progress(console=Console()) as p:
        task = p.add_task("[cyan]Normalizing files...", total=len(files))
        with cf.ThreadPoolExecutor(max_workers=16) as pool:
            futures = {}
            for file in files:
                if file.endswith('.jdx'):
                    file_path = os.path.join(folder_path, file)
                    name = file.split('.')[0]
                    try:
                        futures[pool.submit(normalize_file, file_path, groups[name])] = file
                        p.advance(task)
                    except KeyError:
                        continue
            for future in cf.as_completed(futures):
                file = futures[future]
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error processing {file}: {e}")
            logger.info("All files processed. Shutting down thread pool.")
            pool.shutdown(wait=True)
    print("Normalization complete.")

if __name__ == "__main__":
    folder_path = "/Users/Bram/Downloads/Dataset/nist_IR/IR"
    groups_file = "/Users/Bram/Downloads/Dataset/groups.json"
    logger.info(f"Loading groups from {groups_file}")
    with open(groups_file) as f:
        groups = json.load(f)
    new_groups = dict()
    max_lines = 1
    for group in groups:
        webbook_id = group["filename"].split(".")[0]
        new_groups[webbook_id] = group["groups"]
    logger.info(f"Loaded {len(new_groups)} groups.")
    logger.info(f"Starting normalization in folder: {folder_path}")
    normalize_folder(folder_path, new_groups)