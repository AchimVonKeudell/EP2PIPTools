# -*- coding: utf-8 -*-
"""

@author: Steffen Schüttler
"""

import os
import openpiv as op 
from openpiv import tools, pyprocess, validation, filters, scaling
import numpy as np
import matplotlib.pyplot as plt
import imageio
import pathlib
from typing import Any, Union, List, Optional

# ======= Functions ========

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

def calc_PIV(frame_a,frame_b):
    winsize = 32 # pixels, interrogation window size in frame A
    searchsize = 34 # pixels, search area size in frame B
    overlap = 26 # pixels, 50% overlap
    dt = 0.02 # sec, time interval between the two frames
    
    u0, v0, sig2noise = pyprocess.extended_search_area_piv(
        frame_a.astype(np.int32),
        frame_b.astype(np.int32),
        window_size=winsize,
        overlap=overlap,
        dt=dt,
        search_area_size=searchsize,
        sig2noise_method='peak2peak',
    )
    
    x, y = pyprocess.get_coordinates(
        image_size=frame_a.shape,
        search_area_size=searchsize,
        overlap=overlap,
    )
    
    invalid_mask = op.validation.sig2noise_val(
        u0, v0, sig2noise,
        threshold = 1.4,
    )
    
    u2, v2 = op.filters.replace_outliers(
        u0, v0,
        invalid_mask,
        method='localmean',
        max_iter=5,
        kernel_size=3,
    )
    
    # convert x,y to mm
    # convert u,v to mm/sec
    x, y, u3, v3 = scaling.uniform(
        x, y, u2, v2,
        scaling_factor = 30,  # pixels/millimeter
    )
    
    # 0,0 shall be bottom left, positive rotation rate is counterclockwise
    x, y, u3, v3 = tools.transform_coordinates(x, y, u3, v3)
    
    v_abs = np.sqrt(u3**2+v3**2)
    vector_field = [x, y, u3, v3, invalid_mask[2]]
    
    return vector_field, v_abs 

def plot_PIV(vector_field, v_abs, save_path):

    fig, ax = plt.subplots(figsize=(3.2,6))
    im = ax.imshow(v_abs,cmap="viridis",
                   extent=(np.min(vector_field[0]),np.max(vector_field[0]),np.min(vector_field[1]),np.max(vector_field[1])),interpolation= "bilinear",
                       vmin=0,vmax=3) #,alpha=0.5) #x, y,
    cbar = fig.colorbar(im, ax=ax,shrink=0.8,label="velocity magnitude / mm$\,$s$^{-1}$")
    cbar.ax.tick_params(labelsize=10)
    plt.tight_layout()
    
    fig, ax = display_vector_field_custom( 
        vector_field,
        ax=ax, scaling_factor=30,
        scale=10, # scale defines here the arrow length
        width=0.005, # width is the thickness of the arrow
        on_img=False, # overlay on the image
    ); 
    plt.tight_layout()
    plt.savefig("../figures/" + save_path + ".png")
    # plt.savefig("../figures/" + save_path + ".pdf")
    plt.show()
    
# ======= Select data and execute ========

folder_data = "../data/"
folder_arr = os.listdir(folder_data)
folder_arr = [folder for folder in folder_arr if ".txt" not in folder]
path_dict = {}

for folder in folder_arr:
    file_arr = []
    filenames = os.listdir(folder_data + "/" + folder)
    for file in filenames:
        if file.endswith(".tif"):
            file_arr.append(file)
    path_dict[folder] = file_arr

for key in path_dict.keys():
    print(key)
    if key + "_results" not in os.listdir("../figures/"):
        os.mkdir("../figures/" + key + "_results")

    for number in range(1,6):
        frame_a  = op.tools.imread(folder_data + key + "/" + path_dict[key][int(np.floor(number*len(path_dict[key])/5-2))])
        frame_a = frame_a[0:700,30:290]
        frame_b  = op.tools.imread(folder_data + key + "/" + path_dict[key][int(np.floor(number*len(path_dict[key])/5-1))])
        frame_b = frame_b[0:700,30:290]
        
        PIV, PIV_abs = calc_PIV(frame_a, frame_b)
        plot_PIV(PIV, PIV_abs, key + "_results/" + key + "_" + str(number))