# encoding: utf-8

import os
from .kcorrect import kcorrect, response, template, fitter, utils

__all__ = ["kcorrect", "response", "template", "fitter", "utils"]

NAME = 'kcorrect'

KCORRECT_DIR = os.path.dirname(__file__)

__version__ = '5.1.0b'

