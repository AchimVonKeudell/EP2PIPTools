from functions_PIV_image import display_vector_field_custom,plot_PIV_surface

PIV_measurements = ["PIV_0,5slm","PIV_1,0slm","PIV_2,0slm"]

fig, axes = plt.subplots(figsize=(10,4),nrows=1,ncols=3,constrained_layout=True)
plt.setp(axes, xticks=[-5,-2.5,0,2.5,5],yticks=[0,2.5,5,7.5,10,12.5,15])
axes[0].grid(False)

i = 0
for file_PIV in PIV_measurements:
    print(file_PIV)
    flow = float(re.findall('PIV_(.*)slm', file_PIV)[0].replace(",","."))
    
    if flow == 5:
        flow = flow/10
    
    for file in os.listdir(folder+file_PIV+"/"):
        if re.search('_x', file):
            x = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_y', file):
            y = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_ux', file):
            ux = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_uy', file):
            uy = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_mask', file):
            mask = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_u', file):
            PIV_abs = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
            
    PIV = [x,y,ux,uy,mask]
    PIV[0] = PIV[0]-PIV[0].mean()

    im = axes[i].imshow(PIV_abs,cmap="plasma",
                   extent=(np.max(PIV[0]),np.min(PIV[0]),np.max(PIV[1]),np.min(PIV[1])),interpolation= "bilinear",
                       vmin=0,vmax=3) #,alpha=0.5) #x, y,
    
    display_vector_field_custom( 
        PIV,
        ax=axes[i], scaling_factor=1,
        scale=12, # scale defines here the arrow length
        width=0.005, # width is the thickness of the arrow
        on_img=False, # overlay on the image
    )
    
    axes[i].set_xlim(-5,5)
    axes[i].set_xlabel("$\Delta$x / mm")
    axes[i].set_ylim(np.max(PIV[1]),0)
    axes[i].set_title("%.1f$\,$slm" %flow,fontweight="bold")
    axes[i].grid(False)
    
    i += 1

axes[0].set_ylabel("z / mm")
axes[1].set_ylabel("")
axes[2].set_ylabel("")
cbar = fig.colorbar(im, ax=axes[-1],shrink=1,label="velocity magnitude / mm$\,$s$^{-1}$", pad=0.03)
cbar.ax.tick_params(labelsize=14,size=0, which="both")
# plt.tight_layout()
if True:
    plt.savefig("Results_liquid_phase_PIV.png")
    plt.savefig("Results_liquid_phase_PIV.pdf")
plt.show() 

#%% PIV Zoom

PIV_measurements = ["PIV_0,5slm","PIV_1,0slm","PIV_2,0slm"]

fig, axes = plt.subplots(figsize=(10,2.2),dpi=600,nrows=1,ncols=3,constrained_layout=True)
i = 0
for file_PIV in PIV_measurements:
    print(file_PIV)
    flow = float(re.findall('PIV_(.*)slm', file_PIV)[0].replace(",","."))
    
    if flow == 5:
        flow = flow/10
    
    for file in os.listdir(folder+file_PIV+"/"):
        if re.search('_x', file):
            x = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_y', file):
            y = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_ux', file):
            ux = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_uy', file):
            uy = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_mask', file):
            mask = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
        elif re.search('_u', file):
            PIV_abs = np.array(pd.read_csv(folder+file_PIV+"/"+file,delimiter="\t")).T[1:].T
            
    PIV = [x,y,ux,uy,mask]
    PIV[0] = PIV[0]-PIV[0].mean()

    im = axes[i].imshow(PIV_abs,cmap="plasma",
                   extent=(np.max(PIV[0]),np.min(PIV[0]),np.max(PIV[1]),np.min(PIV[1])),interpolation= "bilinear",
                       vmin=0,vmax=3) #,alpha=0.5) #x, y,
    
    display_vector_field_custom( 
        PIV,
        ax=axes[i], scaling_factor=1,
        scale=12, # scale defines here the arrow length
        width=0.005, # width is the thickness of the arrow
        on_img=False, # overlay on the image
    )
    
    axes[i].set_xlim(-5,5)
    axes[i].set_xlabel("$\Delta$x / mm")
    axes[i].set_ylim(4,0)
    axes[i].set_title("%.1f$\,$slm" %flow,fontweight="bold")
    axes[i].grid(False)
    
    i += 1


axes[0].set_ylabel("z / mm")
axes[1].set_ylabel("")
axes[2].set_ylabel("")

# cbar = fig.colorbar(im, ax=axes,shrink=1,label="velocity magnitude / mm$\,$s$^{-1}$")
# cbar.ax.tick_params(labelsize=14,size=0, which="both")
# plt.tight_layout()
if True:
    plt.savefig("Results_liquid_phase_PIV_zoom.png")
    plt.savefig("Results_liquid_phase_PIV_zoom.pdf")
plt.show()
