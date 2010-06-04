# Copyright (C) 2008 NSC Jyvaskyla
# Please see the accompanying LICENSE file for further information.

from aseinterface import Hotbit
from atoms import Atoms
from atoms import ExtendedTrajectory
from electrostatics import Electrostatics
from element import Element
from elements import Elements
from environment import Environment
from environment import LinearlyPolarizedLaser
from grids import Grids
from interactions import Interactions
from occupations import Occupations
from output import Output
from repulsion import RepulsivePotential
from repulsion import Repulsion
from solver import Solver
from states import States
from wfpropagation import WFPropagation

from hotbit.analysis import LinearResponse
from hotbit.analysis import MullikenAnalysis
from hotbit.analysis import MullikenBondAnalysis
from hotbit.analysis import DensityOfStates 

from hotbit.containers import DoubleChiral
from hotbit.containers import Bravais
from hotbit.containers import Chiral
from hotbit.containers import Wedge
from hotbit.containers import Sphere
from hotbit.containers import Saddle
from hotbit.containers import Gaussian
from hotbit.containers import Slab
from hotbit.containers import ContainerTest1

from hotbit.parametrization import SlaterKosterTable
from hotbit.parametrization import KSAllElectron
from hotbit.parametrization import RepulsiveFitting
from hotbit.parametrization import ParametrizationTest

from box.wedgeatoms import *
from os import environ, path

import atexit
import _hotbit

hbpar = environ.get('HOTBIT_PARAMETERS')
fixpar = path.join(environ.get('HOTBIT_PARAMETERS'),'fixed_parameters')
testpar = path.join(environ.get('HOTBIT_PARAMETERS'),'inofficial')


#
# Free eigenvalue solver workspace on exit
#
atexit.register(_hotbit.free_geig_workspace)
from hotbit.version import hotbit_version

