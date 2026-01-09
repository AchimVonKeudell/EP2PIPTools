"""
This file defines a cache decorator that can be used to cache output of functions, where the function either has no
parameters,a global, or a global and local quanta as parameter. Once the output for this function (for the given quanta)
has been calculated, it will be returned from the cache on any subsequent calls.

Once the reset function is called, the cache is emptied, and will be refilled upon calling the functions.

NOTE: Don't forget to empty the cache when caching distribution function calculations once the temperature of the
distribution is changed.
"""
import functools


def reset(original_function=None, check_parameter_change=False):
    """
    Reset the cache upon calling the given function
    :return:
    """
    def _decorate(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            function_id = id(func)
            function_vars = vars()

            if not hasattr(self, '_cache_reset_parameters'):
                self._cache_reset_parameters = {}

            prev_func_params = self._cache_reset_parameters.setdefault(function_id, {})
            if not check_parameter_change or prev_func_params != function_vars:
                if hasattr(self, '_cache'):
                    self._cache.setdefault('global_quanta', {}).clear()
                    self._cache.setdefault('fixed', {}).clear()

            self._cache_reset_parameters[function_id] = function_vars

            func(self, *args, **kwargs)

        return wrapper

    if original_function:
        return _decorate(original_function)

    return _decorate


def cache_by_global_quanta_decorator(original_function=None, immutable=False):
    """
    Cache the output of a function that is calculated for a given global quanta
    :param func:
    :return:
    """

    def _decorate(func):
        @functools.wraps(func)
        def wrapper(self, global_quanta):
            function_id = id(func)

            if immutable:
                key = 'global_quanta_immutable'
            else:
                key = 'global_quanta'

            try:
                return self._cache[key][function_id][global_quanta.hitran_raw_quanta]
            except (AttributeError, KeyError):
                if not hasattr(self, '_cache'):
                    self._cache = {}

                self._cache.setdefault(key, {})
                self._cache[key].setdefault(function_id, {})
                self._cache[key][function_id][global_quanta.hitran_raw_quanta] = func(self, global_quanta)
                return self._cache[key][function_id][global_quanta.hitran_raw_quanta]

        return wrapper

    if original_function:
        return _decorate(original_function)

    return _decorate


def cache_fixed_value(original_function=None, immutable=False):
    """
    Cache the output of a function that is calculated for a given global quanta
    :param func:
    :return:
    """
    def _decorate(func):
        @functools.wraps(func)
        def wrapper(self):
            function_id = id(func)

            if immutable:
                key = 'fixed_immutable'
            else:
                key = 'fixed'

            try:
                return self._cache[key][function_id]
            except (AttributeError, KeyError):
                if not hasattr(self, '_cache'):
                    self._cache = {}

                self._cache.setdefault(key, {})
                self._cache[key].setdefault(function_id, {})
                self._cache[key][function_id] = func(self)
                return self._cache[key][function_id]

        return wrapper

    if original_function:
        return _decorate(original_function)

    return _decorate
