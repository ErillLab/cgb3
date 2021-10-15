"""The module for miscellaneous functions."""

import os
import math
import numpy as np
import tempfile


def mean(xs):
    """Returns the average of a given list of numbers."""
    return sum(xs) / float(len(xs))


def log2(x):
    """Returns the logarithm in base 2."""
    return math.log(x, 2)


def normalize(xs):
    """Normalize the given list of values."""
    return [float(x)/sum(xs) for x in xs]


def weighted_choice(xs, weights, count=1):
    """Given a list of elements, randomly selects one based on weights.
    Args:
        xs (list): elements
        weights (list): weights
    """
    ps = [float(w)/sum(weights) for w in weights]
    return np.random.choice(xs, size=count, p=ps)


def unique(xs, f):
    """Makes the list xs unique by comparing elements with f(x).

    Keeps the first occurrence of each element.
    """
    unique_list = []
    for x in xs:
        if not f(x) in list(map(f, unique_list)):
            unique_list.append(x)
    return unique_list


def temp_file_name(dir=tempfile.gettempdir(), prefix='', suffix=''):
    """Creates a unique file and return its name."""
    tmpfile = tempfile.NamedTemporaryFile(dir=dir, prefix=prefix,
                                          suffix=suffix, delete=False)
    name = tmpfile.name
    tmpfile.close()
    return name


def directory(*paths):
    """Makes a directory specified by one or more args.

    Creates the directory if it doesn't exist.

    Args:
        paths ([string]): the full path for the directory
    Returns:
        string: the full path for the directory
    """
    d = os.path.join(*paths)
    if not os.path.exists(d):
        os.makedirs(d)
    return d
