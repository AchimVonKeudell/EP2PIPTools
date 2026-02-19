# -*- coding: utf-8 -*-
"""
Simulation von Diffusion + Konvektion von H2O2
mit Abschalten des Quellterms nach bestimmter Zeit.

Autor: Steffen Schüttler, Sept 2025
"""

import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import shutil
from tqdm import tqdm

plt.rcParams.update({'font.size': 12})

# =============================================================================
# Globale Parameter
# =============================================================================
ZNB = 220          # Anzahl Gitterpunkte
GP = 1             # Ghost Points
DZ = 1e-4          # Zellgröße [m]
N0 = 2.4e28        # Atomdichte [1/m^3]

DT0 = 5e-4         # Startzeit-Schritt [s]
T_END = 1800       # Endzeit [s]
DT_PLOT = 1        # Zeitschritt für Speicherung/Plot [s]

J_H2O2 = 2.1e20    # Fluss [1/m^2/s]
DIFF_COEFF = 2.6e-8  # Diffusionskoeffizient [m^2/s]
C_MAX = 3e26       # Max. Konzentration an Oberfläche
U = 3e-3           # Konvektionsgeschwindigkeit [m/s]

SOURCE_OFF = 1200  # Zeitpunkt, ab dem die Quelle abgeschaltet wird [s]


# =============================================================================
# Hilfsfunktionen
# =============================================================================
def initialize_grid(znb: int, gp: int, dz: float) -> np.ndarray:
    return np.linspace(0, (znb + 2 * gp - 1) * dz, znb + 2 * gp)


def apply_boundary_conditions(conc: np.ndarray, gp: int) -> None:
    conc[:gp] = conc[gp]
    conc[-gp:] = conc[-gp - 1]


def compute_derivative(conc: np.ndarray, diff: float, dz: float,
                       flux: float, gp: int, velocity: float,
                       influx_on: bool = True) -> np.ndarray:
    deriv = np.zeros_like(conc)
    for i in range(gp, len(conc) - gp):
        inflow = flux / dz if (influx_on and i == 1) else 0.0
        diffusion = diff * ((conc[i + 1] - conc[i]) +
                            (conc[i - 1] - conc[i])) / dz**2
        convection = -velocity * (conc[i + 1] - conc[i - 1]) / (2 * dz)
        deriv[i] = diffusion + convection + inflow
    return deriv


def euler_step(conc: np.ndarray, deriv: np.ndarray,
               dt: float, dz: float, diff: float) -> None:
    dt_diff = 0.5 * dz**2 / diff
    dt_eff = min(dt, dt_diff)
    conc += deriv * dt_eff
    conc[conc < 0] = 0


# =============================================================================
# Simulation
# =============================================================================
def run_simulation(diff: float, velocity: float,
                   save_history: bool = True) -> dict:
    z = initialize_grid(ZNB, GP, DZ)
    conc = np.zeros_like(z)
    conc[1] = 1.0

    times, averages, profiles = [], [], []

    n_steps = int(T_END / DT_PLOT)
    for step in tqdm(range(n_steps), desc=f"Simulation (u={velocity:.2e} m/s)"):
        t = 0.0
        while t < DT_PLOT:
            influx_on = (step * DT_PLOT) < SOURCE_OFF
            deriv = compute_derivative(conc, diff, DZ, J_H2O2, GP,
                                       velocity, influx_on)
            euler_step(conc, deriv, DT0, DZ, diff)
            apply_boundary_conditions(conc, GP)
            t += DT0

        times.append(step * DT_PLOT)
        averages.append(np.mean(conc))
        if save_history:
            profiles.append(conc.copy())

    return {
        "depth": z,
        "conc": conc,
        "time": np.array(times),
        "averages": np.array(averages),
        "profiles": np.array(profiles) if save_history else None
    }


# =============================================================================
# Animation
# =============================================================================
def animate_profiles(results, save_path="simulation.mp4", fps=15):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    line1, = ax1.plot([], [], lw=2, color="blue")
    line2, = ax2.plot([], [], lw=2, color="black")

    ax1.set_xlim(0, results["depth"].max() * 1e3)
    ax1.set_ylim(0, results["profiles"].max() * 1.1)
    ax1.set_xlabel("Depth (mm)")
    ax1.set_ylabel("Concentration (m$^{-3}$)")

    ax2.set_xlim(0, results["time"][-1])
    ax2.set_ylim(0, max(results["averages"]) * 1.1)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Average concentration")

    ax2.axvline(x=SOURCE_OFF, color="red", linestyle="--",
                label=f"Source off ({SOURCE_OFF:.0f} s)")
    ax2.legend()

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        return line1, line2

    def update(frame):
        line1.set_data(results["depth"] * 1e3, results["profiles"][frame])
        line2.set_data(results["time"][:frame + 1],
                       results["averages"][:frame + 1])
        fig.suptitle(f"t = {results['time'][frame]:.1f} s")
        return line1, line2

    ani = animation.FuncAnimation(
        fig, update, frames=len(results["time"]),
        init_func=init, blit=True, interval=1000 / fps
    )

    if save_path.endswith(".mp4"):
        if shutil.which("ffmpeg") is None:
            raise RuntimeError("Fehler: ffmpeg nicht gefunden!")
        writer = animation.FFMpegWriter(fps=fps)
    elif save_path.endswith(".gif"):
        if shutil.which("magick") is None and shutil.which("convert") is None:
            raise RuntimeError("Fehler: ImageMagick nicht gefunden!")
        writer = animation.ImageMagickWriter(fps=fps)
    else:
        raise ValueError("Ungültige Endung: Bitte .mp4 oder .gif angeben")

    ani.save(save_path, writer=writer, dpi=150)
    plt.close(fig)
    print(f"✅ Animation gespeichert: {save_path}")


# =============================================================================
# Vergleichsdiagramm mit Heatmaps
# =============================================================================
def plot_comparison(results_diff, results_conv, save_path="comparison.png"):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax1, ax2, ax3, ax4 = axes.flatten()

    # --- Konzentrationsprofile ---
    times_to_plot = [600, 1200, 1800]
    colors = ["blue", "green", "purple"]

    for t_plot, color in zip(times_to_plot, colors):
        idx_d = np.argmin(np.abs(results_diff["time"] - t_plot))
        idx_c = np.argmin(np.abs(results_conv["time"] - t_plot))

        ax1.plot(results_diff["depth"] * 1e3, results_diff["profiles"][idx_d],
                 label=f"Diffusion {t_plot}s", color=color, linestyle="--")
        ax1.plot(results_conv["depth"] * 1e3, results_conv["profiles"][idx_c],
                 label=f"Conv.+Diff. {t_plot}s", color=color)

    ax1.set_xlabel("Depth (mm)")
    ax1.set_ylabel("Concentration (m$^{-3}$)")
    ax1.legend()
    ax1.set_title("Konzentrationsprofile")

    # --- Mittelwertkurven ---
    ax2.plot(results_diff["time"], results_diff["averages"],
             label="Diffusion", color="blue", linestyle="--")
    ax2.plot(results_conv["time"], results_conv["averages"],
             label="Conv.+Diff.", color="black")

    ax2.axvline(x=SOURCE_OFF, color="red", linestyle="--",
                label=f"Source off ({SOURCE_OFF:.0f} s)")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Average concentration")
    ax2.legend()
    ax2.set_title("Mittlere Konzentration")

    # --- Heatmap: Diffusion ---
    im3 = ax3.imshow(results_diff["profiles"].T,
                     aspect="auto", origin="lower",
                     extent=[results_diff["time"][0], results_diff["time"][-1],
                             results_diff["depth"].min() * 1e3,
                             results_diff["depth"].max() * 1e3],
                     cmap="viridis")
    ax3.set_xlabel("Time (s)")
    ax3.set_ylabel("Depth (mm)")
    ax3.set_title("Heatmap Diffusion")
    fig.colorbar(im3, ax=ax3, label="Concentration")

    # --- Heatmap: Konvektion + Diffusion ---
    im4 = ax4.imshow(results_conv["profiles"].T,
                     aspect="auto", origin="lower",
                     extent=[results_conv["time"][0], results_conv["time"][-1],
                             results_conv["depth"].min() * 1e3,
                             results_conv["depth"].max() * 1e3],
                     cmap="viridis")
    ax4.set_xlabel("Time (s)")
    ax4.set_ylabel("Depth (mm)")
    ax4.set_title("Heatmap Conv.+Diff.")
    fig.colorbar(im4, ax=ax4, label="Concentration")

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"✅ Vergleichsdiagramm mit Heatmaps gespeichert: {save_path}")


# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    print("Simulation gestartet:", datetime.datetime.now())

    results_diff = run_simulation(DIFF_COEFF, velocity=0.0, save_history=True)
    results_conv = run_simulation(DIFF_COEFF, velocity=U, save_history=True)

    animate_profiles(results_diff, save_path="diffusion.mp4", fps=15)    
    animate_profiles(results_conv, save_path="convection_diffusion.mp4", fps=15)
    plot_comparison(results_diff, results_conv, save_path="comparison.png")

    print("Fertig:", datetime.datetime.now())
