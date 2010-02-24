# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from aseinterface import Calculator, Hotbit
from element import Element
from repulsion import RepulsivePotential
from wfpropagation import WFPropagation
from environment import Environment, LinearlyPolarizedLaser
from output import Output
from states import States
from box.wedgeatoms import *
from hotbit.atoms import Atoms
from os import environ, path

from hotbit.parametrization.slako import SlaterKosterTable
from hotbit.parametrization.atom import KSAllElectron
from hotbit.parametrization.fitting import RepulsiveFitting

fixpar = path.join(environ.get('HOTBIT_PARAMETERS'),'fixed_parameters')
testpar = path.join(environ.get('HOTBIT_PARAMETERS'),'inofficial')
