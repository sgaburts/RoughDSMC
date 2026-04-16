# Renames SPARTA variable in .vtu and .pvtu files in a directory
import xml.etree.ElementTree as ET
import sys
import os
import glob

# Define grouped renaming
grouped_rename_map = {
    ("f_fgrid[1]", "f_fgrid2[1]", "f_fgrid3[1]"): "nrho",
    ("f_fgrid[2]", "f_fgrid2[2]", "f_fgrid3[2]"): "u",
    ("f_fgrid[3]", "f_fgrid2[3]", "f_fgrid3[3]"): "v",
    ("f_fgrid[4]", "f_fgrid2[4]", "f_fgrid3[4]"): "w",
    ("f_fgrid[5]", "f_fgrid2[5]", "f_fgrid3[5]"): "temp",
    ("f_fgrid[6]", "f_fgrid2[6]", "f_fgrid3[6]"): "T_vib",
    ("f_fgrid[7]", "f_fgrid2[7]", "f_fgrid3[7]"): "T_rot",
    ("f_fgrid[8]", "f_fgrid2[8]", "f_fgrid3[8]"): "n",
    ("f_fgrid[9]", "f_fgrid2[9]", "f_fgrid3[9]"): "massrho",
    ("f_fgrid[10]", "f_fgrid2[10]", "f_fgrid3[10]"): "press",
    ("f_fgrid[11]", "f_fgrid2[11]", "f_fgrid3[11]"): "ke",
    ("f_fgrid[12]", "f_fgrid2[12]", "f_fgrid3[12]"): "kerho",
    ("c_cgrid[7]", "c_cgrid2[7]", "c_cgrid3[7]"): "n_inst",
    ("f_forces_circle2[1]","f_forces_circle2[1]","f_forces_circle2[1]"): "n",
    ("f_forces_circle2[2]","f_forces_circle2[2]","f_forces_circle2[2]"): "press",
    ("f_forces_circle2[3]","f_forces_circle2[3]","f_forces_circle2[3]"): "fx",
    ("f_forces_circle2[4]","f_forces_circle2[4]","f_forces_circle2[4]"): "fy",
    ("f_forces_circle2[5]","f_forces_circle2[5]","f_forces_circle2[5]"): "fz",
    ("f_forces_circle2[6]","f_forces_circle2[6]","f_forces_circle2[6]"): "px",
    ("f_forces_circle2[7]","f_forces_circle2[7]","f_forces_circle2[7]"): "py",
    ("f_forces_circle2[8]","f_forces_circle2[8]","f_forces_circle2[8]"): "pz",
    ("f_forces_circle2[9]","f_forces_circle2[9]","f_forces_circle2[9]"): "shx",
    ("f_forces_circle2[10]","f_forces_circle2[10]","f_forces_circle2[10]"): "shy",
    ("f_forces_circle2[11]","f_forces_circle2[11]","f_forces_circle2[11]"): "shz",
    ("f_forces_circle2[12]","f_forces_circle2[12]","f_forces_circle2[12]"): "mflux",
    ("c_forces_circle[1]","c_forces_circle[1]","c_forces_circle[11]"): "n_inst",
}


rename_map = {old: new for group, new in grouped_rename_map.items() for old in group}

def rename_dataarrays(tree, tagname):
    modified = False
    for data_array in tree.getroot().iter(tagname):
        name = data_array.get("Name")
        if name in rename_map:
            data_array.set("Name", rename_map[name])
            print(f"  Renamed '{name}' → '{rename_map[name]}'")
            modified = True
    return modified

def process_vtu_file(filepath):
    try:
        tree = ET.parse(filepath)
        if rename_dataarrays(tree, tagname="DataArray"):
            tree.write(filepath, encoding="utf-8", xml_declaration=True)
            print(f"Updated .vtu: {filepath}")
        else:
            print(f"No changes in .vtu: {filepath}")
    except Exception as e:
        print(f"Failed to process {filepath}: {e}")

def process_pvtu_file(filepath):
    try:
        tree = ET.parse(filepath)
        root = tree.getroot()

        # Extract referenced .vtu files first
        vtu_relative_files = [e.get("Source") for e in root.iter("Piece") if e.get("Source")]
        pvtu_dir = os.path.dirname(os.path.abspath(filepath))
        vtu_files = [os.path.join(pvtu_dir, f) for f in vtu_relative_files]

        # Rename PDataArray entries
        if rename_dataarrays(tree, tagname="PDataArray"):
            tree.write(filepath, encoding="utf-8", xml_declaration=True)
            print(f"Updated .pvtu: {filepath}")
        else:
            print(f"No changes in .pvtu: {filepath}")

        # Process referenced .vtu files
        for vtu in vtu_files:
            print(f"Processing referenced .vtu: {vtu}")
            if os.path.isfile(vtu):
                process_vtu_file(vtu)
            else:
                print(f"File not found: {vtu}")
    except Exception as e:
        print(f"Failed to process {filepath}: {e}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python rename_variables.py <directory_or_file>")
        sys.exit(1)

    path = sys.argv[1]

    if os.path.isdir(path):
        # Process all .pvtu and .vtu files in directory
        for filepath in glob.glob(os.path.join(path, "*.pvtu")):
            print(f"\nFound .pvtu: {filepath}")
            process_pvtu_file(filepath)
        for filepath in glob.glob(os.path.join(path, "*.vtu")):
            print(f"\nFound .vtu: {filepath}")
            process_vtu_file(filepath)
    elif os.path.isfile(path):
        if path.endswith(".pvtu"):
            process_pvtu_file(path)
        elif path.endswith(".vtu"):
            process_vtu_file(path)
        else:
            print(f"Unsupported file type: {path}")
    else:
        print(f"Invalid path: {path}")

if __name__ == "__main__":
    main()
