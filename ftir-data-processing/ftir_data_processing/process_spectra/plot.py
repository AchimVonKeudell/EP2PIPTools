import matplotlib.pyplot as plt


def plot_transmittance_spectrum(filepath_figure, wavenumbers, original_spectrum, corrected_spectrum, **kwargs):
    """Plotting a transmittance spectrum

    :param filepath_figure:
    :param wavenumbers:
    :param original_spectrum:
    :param corrected_spectrum:

    """
    if 'save_space' in kwargs:
        if kwargs['save_space']:
            return

    if 'plot_data' in kwargs:
        if not kwargs['plot_data']:
            return

    plt.rc('axes', labelweight='bold')

    fig = plt.figure(figsize=(7, 4))
    ax = fig.add_subplot()
    ax.set(xlabel='wavenumbers [cm$^{-1}$]', xlim=(wavenumbers[0], wavenumbers[-1]),
           ylabel=r'transmittance [-]')
    if 'y_plot_limits' in kwargs:
        ax.set(ylim=kwargs['y_plot_limits'])
    ax.grid(True)

    ax.plot(wavenumbers, original_spectrum, label='measured spectrum')
    if 'baseline' in kwargs:
        ax.plot(wavenumbers, kwargs['baseline'], label='baseline', c='black', ls='--')
    if 'best_fit_air_contributions' in kwargs:
        best_fit_air_contributions = kwargs['best_fit_air_contributions']
        [ax.plot(wavenumbers, y, label=label) for label, y in best_fit_air_contributions.items()]
    ax.plot(wavenumbers, corrected_spectrum, label='corrected spectrum', c='black')

    ax.legend()
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.01)

    if filepath_figure != '':
        fig.savefig(filepath_figure)
        print("Saved:", filepath_figure)
        plt.close()
    return