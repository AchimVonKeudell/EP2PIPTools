# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:24:27 2023

@author: Steffen
"""

import numpy as np
import matplotlib.pyplot as plt

from openpiv import tools
import imageio
from typing import Any, Union, Optional
import pathlib


def display_vector_field_custom(
    vector_field, #filename: Union[pathlib.Path, str],
    on_img: Optional[bool]=False,
    image_name: Optional[Union[pathlib.Path,str]]=None,
    window_size: Optional[int]=32,
    scaling_factor: Optional[float]=1.,
    ax: Optional[Any]=None,
    width: Optional[float]=0.0025,
    show_invalid: Optional[bool]=True,
    **kw
):

    # print(f' Loading {filename} which exists {filename.exists()}')
    a = vector_field
    # parse
    x, y, u, v, mask = a[0], a[1], a[2], a[3], a[4] #, a[:, 5]


    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    if on_img is True:  # plot a background image
        im = imageio.imread(image_name)
        im = tools.negative(im)  # plot negative of the image for more clarity
        xmax = np.amax(x) + window_size / (2 * scaling_factor)
        ymax = np.amax(y) + window_size / (2 * scaling_factor)
        ax.imshow(im, cmap="Greys_r", extent=[0.0, xmax, 0.0, ymax])    


    # first mask whatever has to be masked
    u[mask.astype(bool)] = 0.
    v[mask.astype(bool)] = 0.
    
    # now mark the valid/invalid vectors
    invalid = mask > 0 
    valid = ~invalid

    ax.quiver(
        x[valid],
        y[valid],
        u[valid],
        v[valid],
        color="w",
        width=width,
        **kw
        )
        
    if show_invalid and len(invalid) > 0:
        ax.quiver(
                x[invalid],
                y[invalid],
                u[invalid],
                v[invalid],
                color="r",
                width=width,
                **kw,
                )
    ax.set_aspect(1.)
    # fig.canvas.set_window_title('Vector field, '+str(np.count_nonzero(invalid))+' wrong vectors')

    ax.set_xlabel("width / mm")
    ax.set_ylabel("height / mm")
    # ax.set_xlim(0)
    # ax.set_ylim(0)

    return fig, ax

def plot_PIV_surface(vector_field, v_abs, colormap, flow, save_path, SAVE):

    fig, ax = plt.subplots(figsize=(3.2,3), dpi=300)
    im = ax.imshow(v_abs,cmap=colormap,
                   extent=(np.max(vector_field[0]),np.min(vector_field[0]),np.max(vector_field[1]),np.min(vector_field[1])),interpolation= "bilinear",
                       vmin=0,vmax=3) #,alpha=0.5) #x, y,
    cbar = fig.colorbar(im, ax=ax,shrink=1,label="velocity magnitude / mm$\,$s$^{-1}$")
    cbar.ax.tick_params(labelsize=10)

    fig, ax = display_vector_field_custom( 
        vector_field,
        ax=ax, scaling_factor=1,
        scale=12, # scale defines here the arrow length
        width=0.005, # width is the thickness of the arrow
        on_img=False, # overlay on the image
    );
    
    ax.set_xlim(0,10)
    # ax.invert_yaxis()
    ax.set_ylim(np.max(vector_field[1]),0)
    # plt.title("v_abs:\n" + "mean: " + str(np.mean(v_abs).round(2)) + "; max: " + str(np.max(v_abs).round(2)), size=10)
    plt.title("PIV - %.1f$\,$slm" %flow,fontsize=10)
    plt.tight_layout()
    if SAVE:
        plt.savefig("../figures/" + save_path + ".png",dpi=300)
        plt.savefig("../figures/" + save_path + ".pdf",dpi=300)
    plt.show()
