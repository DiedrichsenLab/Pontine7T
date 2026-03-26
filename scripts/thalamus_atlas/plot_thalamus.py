import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from nilearn.plotting import plot_roi
from nibabel import Nifti1Image
import matplotlib.patches as mpatches
from distinctipy import get_colors
from distinctipy import distinctipy


def plot_thalamus_wta(
        wta_data,           # 3D array of integer labels
        bg_img=None,        # Optional Nifti background
        mask_img= None,
        lut_cmap=None,     # Nx3 array of colors from LUT
        names=None,
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

    if bg_img is not None:
        affine = bg_img.affine
    elif mask_img is not None:
        affine = mask_img.affine
    else:
        affine = np.eye(4)

    # Ensure colormap is a ListedColormap
    if lut_cmap is None:
        raise ValueError("You must provide a LUT colormap array (Nx3)")
    
    colors = get_colors(len(names))
    color_dict = {name: colors[i] for i, name in enumerate(names)}
    cmap_colors = np.array([color_dict[name] for name in names])
    cmap = ListedColormap(cmap_colors)

    
    #colors = get_colors(len(names), seed=42)  # generates n visually distinct RGB tuples
    #color_dict = {name: colors[i] for i, name in enumerate(names)}
    #cmap_colors = np.array([color_dict[name] for name in names])
    #cmap = ListedColormap(np.array(colors))

    #cmap = ListedColormap(lut_cmap / 255)  # normalize RGB 0-1

    # Convert WTA array to Nifti image
    if bg_img is not None:
        affine = bg_img.affine
    else:
        affine = mask_img.affine

    wta_img = Nifti1Image(wta_data.astype(np.int32), affine)

    # Create figure and axes
    plt.close('all')
    n_slices = len(z_coords)
    fig, axes = plt.subplots(n_slices, 1, figsize=figsize)
    axes = np.atleast_1d(axes)
    plt.figure(fig.number)

    for ax, z in zip(axes, z_coords):
        display = plot_roi(
        wta_img,
        bg_img=bg_img,
        cmap=cmap,
        display_mode='z',
        cut_coords=[z],
        axes=ax,         
        figure=fig,       
        black_bg=True,
        annotate=False
    )
    patches = [mpatches.Patch(color=cmap.colors[i], label=names[i]) for i in range(len(names))]
    fig.legend(handles=patches, loc='right', bbox_to_anchor=(1.0, 0.5), borderaxespad=0., fontsize=8)
    fig.subplots_adjust(right=0.75)

#plt.tight_layout()
#plt.tight_layout(rect=[0, 0, 0.1, 1])