import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from nilearn.plotting import plot_roi
from nibabel import Nifti1Image

def plot_thalamus_wta(
        wta_data,           # 3D array of integer labels
        bg_img=None,        # Optional Nifti background
        lut_cmap=None,      # Nx3 array of colors from LUT
        z_coords=[-7,-5,-3,1,8,12,15],
        figsize=(6, 10)
    ):
    """
    Plot a categorical winner-takes-all thalamus map using a LUT.

    Args:
        wta_data (ndarray): 3D volume of integer labels (1..K)
        bg_img (Nifti1Image or None): Background image (structural)
        lut_cmap (ndarray): Nx3 array of RGB colors (0-255)
        z_coords (list): list of z slices to plot
        figsize (tuple): figure size

    Returns:
        fig (plt.Figure), axes (list)
    """

    # Ensure colormap is a ListedColormap
    if lut_cmap is None:
        raise ValueError("You must provide a LUT colormap array (Nx3)")

    cmap = ListedColormap(lut_cmap / 255)  # normalize RGB 0-1

    # Convert WTA array to Nifti image
    if bg_img is not None:
        affine = bg_img.affine
    else:
        affine = np.eye(4)

    wta_img = Nifti1Image(wta_data.astype(np.int32), affine)

    # Create figure and axes
    n_slices = len(z_coords)
    fig, axes = plt.subplots(n_slices, 1, figsize=figsize)
    axes = np.atleast_1d(axes)

    # Plot each z-slice
    for i, z in enumerate(z_coords):
        plot_roi(
            wta_img,
            bg_img=bg_img,
            cmap=cmap,
            display_mode='z',
            cut_coords=[z],
            axes=axes[i],
            black_bg=True,
            annotate=False
        )

    plt.tight_layout()
    return fig, axes