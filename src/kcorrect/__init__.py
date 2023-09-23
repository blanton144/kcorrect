# encoding: utf-8

import os

KCORRECT_DIR = os.path.dirname(__file__)

from . import fitter, kcorrect, response, template, utils

__all__ = ["kcorrect", "response", "template", "fitter", "utils"]

NAME = 'kcorrect'


try:
    from kcorrect._version import __version__, __version_tuple__
except:
    __version__ = 'work'
