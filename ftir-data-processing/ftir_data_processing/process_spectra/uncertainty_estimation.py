from numpy import ndarray, sum, square


def uncertainty_analysis(x_values: ndarray, y_values: ndarray, function: callable, best_fit_parameter_value: ndarray,
                         internal_error: float, number_of_parameters: int = 2, increment: float = 0.01, *,
                         required_increase_in_reduced_chi_squared=1.0, verbose=False):
    """Estimates the accuracy of the parameter (defined with parameter_name) using the reduced chi squared test, where
            red_chi_sq(parameter_best_fit) + 1 = red_chi_sq(parameter_best_fit + error_parameter)
     is solved for error_parameter.

    :param x_values:
    :param y_values:
    :param function:    callable in the form of f(x, value), i.e. with two variables
    :param best_fit_parameter_value: value found previously with e.g. scipy.optimise.curve_fit
    :param internal_error:
    :param number_of_parameters: # of fit parameters (i.e. which are optimised)
    :param increment:
    :param required_increase_in_reduced_chi_squared:
    :param verbose: boolean to indicate if results must be printed
    :return;
    """
    if not callable(function):
        raise ValueError('function should be callable')

    normalising_red_chi_sq_factor = 1 / (square(internal_error) * (y_values.size - number_of_parameters))

    def reduced_chi_squared(parameter_value) -> float:
        return normalising_red_chi_sq_factor * sum(square(y_values - function(x_values, parameter_value)))

    def get_best_fit_value():
        if isinstance(best_fit_parameter_value, ndarray):
            return best_fit_parameter_value.copy()
        elif isinstance(best_fit_parameter_value, float):
            return best_fit_parameter_value
        raise ValueError(f'What is the type of {best_fit_parameter_value=}??')

    chi_sq_best_fit = reduced_chi_squared(best_fit_parameter_value)
    if verbose:
        print(f'\tmin(red. chi^2): {chi_sq_best_fit:.2f}')

    error_margins = []
    for sign_of_increment in [-1, 1]:
        parameter_set_value = get_best_fit_value()
        increase_chi_sq = 0.

        while (increase_chi_sq < required_increase_in_reduced_chi_squared) & (parameter_set_value > 0.):
            parameter_set_value += sign_of_increment * increment
            if parameter_set_value < 0:
                parameter_set_value = 0.
            increase_chi_sq = reduced_chi_squared(parameter_set_value) - chi_sq_best_fit

        error_margins.append(parameter_set_value - best_fit_parameter_value)

    return chi_sq_best_fit, tuple(error_margins)


