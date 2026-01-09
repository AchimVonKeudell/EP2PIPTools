import matplotlib.pyplot as plt
import numpy as np


def plot_spectra(save_file, wavenumbers, transmittance, fit, **kwargs):

    plt.rc('axes', labelweight='bold')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 4), gridspec_kw={'height_ratios': [3, 1]}, sharex='col')
    ax1.set(xlim=(wavenumbers.min(), wavenumbers.max()), xlabel='wavenumbers/cm$^{-1}$', ylabel='transmittance / -')

    if 'y_plot_limits' in kwargs:
        ax1.set(ylim=kwargs['y_plot_limits'])
        ax2.set(ylim=[value-1 for value in kwargs['y_plot_limits']])

    ax2.set(ylabel='residual / -')
    ax1.grid(True), ax2.grid(True)

    ax1.plot(wavenumbers, transmittance, label='data')
    ax1.plot(wavenumbers, fit, label='fit')

    residual = transmittance - fit

    if 'baseline_spectrum' in kwargs:
        residual -= kwargs['baseline_spectrum']
        ax1.plot(wavenumbers, kwargs['baseline_spectrum'], ls='--', c='black', label='baseline')

    if 'zoom' in kwargs:
        residual *= kwargs['zoom']
        ax2.set(ylabel=f"residual / {kwargs['zoom']:.1e}")
    ax2.plot(wavenumbers, residual)


    ax1.legend(loc=3)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.01)

    if save_file != '':
        fig.savefig(save_file)
        print("Saved:", save_file)
        plt.close()
    return



