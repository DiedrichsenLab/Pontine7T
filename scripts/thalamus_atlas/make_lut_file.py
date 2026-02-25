import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

wk_dir = '/Users/incehusain/fs_projects'

df = pd.read_csv(f"{wk_dir}/thalamus_indices_47_lut.csv")
names = df["name"].tolist()

K = len(names)
lid = np.arange(K)

# --- Extract base parcel names (remove hemisphere prefix)
def get_base_name(name):
    if name.startswith("Left-"):
        return name.replace("Left-", "")
    elif name.startswith("Right-"):
        return name.replace("Right-", "")
    else:
        return name  # e.g. "0"

base_names = [get_base_name(n) for n in names]

# Unique parcel identities (ignoring hemisphere)
unique_bases = sorted(set(base_names))

# Remove background if present
if 0 in unique_bases:
    unique_bases.remove("0")

# --- Generate distinct colors for unique base parcels
n_unique = len(unique_bases)
cmap = plt.get_cmap("hsv", n_unique)
base_colors = (cmap(np.arange(n_unique))[:, :3] * 255).astype(int)

# Map base name -> color
color_dict = dict(zip(unique_bases, base_colors))

# --- Build full color array
colors = np.zeros((K, 3), dtype=int)

for i, base in enumerate(base_names):
    if base == "0":
        colors[i] = [0, 0, 0]  # background
    else:
        colors[i] = color_dict[base]

# --- Write LUT
lut_path = f"{wk_dir}/thalamus_atlas.lut"
with open(lut_path, "w") as f:
    for i in range(K):
        r, g, b = colors[i]
        f.write(f"{lid[i]} {r} {g} {b} {names[i]}\n")

