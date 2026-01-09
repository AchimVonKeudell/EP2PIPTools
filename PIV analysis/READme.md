This script can be used to analyse images of particle image velocimetry (PIV) measurements. It is based on the openPIV packages available for Python. More information about openPIV for Python can be found in https://openpiv.readthedocs.io/en/latest/src/piv_basics.html.

The script "PIV_evaluation_script.py" analyses the raw camera images and creates tables containing the mask (position of valid analysis points), the velocity magnitude, the x and y components of the velocity field and position x and y. The script "plot_PIV_image.py" uses the tables and creates the final PIV image. A zoom to certain positions in the PIV image is also possible. "functions_PIV_image.py" contains the necessary functions to plot the PIV images.

<img width="6000" height="2400" alt="Results_liquid_phase_PIV" src="https://github.com/user-attachments/assets/5efcf397-9c47-41a9-b735-4014987656c0" />

<img width="6000" height="1320" alt="Results_liquid_phase_PIV_zoom" src="https://github.com/user-attachments/assets/ad133025-a75d-462f-a056-9b9f816a17bc" />

Images taken from Steffen Schüttler, Dissertation, Ruhr University Bochum, 2025
