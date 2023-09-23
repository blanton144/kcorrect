# encoding: utf-8

import os

KCORRECT_DIR = os.path.dirname(__file__)

from . import fitter, kcorrect, response, template, utils

__all__ = ["kcorrect", "response", "template", "fitter", "utils"]

NAME = 'kcorrect'

__version__ = '5.1.0b'

