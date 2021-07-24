from __future__ import absolute_import, division, print_function

from .common.frame import Frame
from .common.block import Block
from .common.distance import Distance
from .common.com_distance import COMDistance
from .common.pbc import PBC
from .common.traj_to_numpy import TrajectoryToNumpy
from .common.units import *
from .common.wholemolecules import Wholemolecules
from .common.cluster import Cluster
from .common.fit import Fit
#from .compute.packing_defects import PackingDefects
from .compute.pd import PackingDefect, PackingDefectPMDA
from .compute.pd2 import PackingDefect2, PackingDefect2PMDA
from .compute.order_parameter import OrderParameters
from .compute.rdf import RDF
from .compute.density_new import Density
#from .compute.density_help import density_help
#from .compute.density import Density
#from .compute.old_density import OldDensity
from .compute.membrane import Membrane
from .compute.coord import Coord
from .compute.dihedral import Dihedral
from .common.covariance import Covariance
#from .scripts.reus_analysis import REUS_Analysis
#from .scripts.lammpsdatawriter import LAMMPSDATAWriter
from .scripts.LAMMPS           import DATAWriter
from .scripts.lammpstrjwriter import LAMMPSTRJWriter
from .scripts.lammpstrjreader import LAMMPSTRJReader
from .scripts.cgmapping import CGMapping, CGMappingPMDA
from .scripts.cgmapping_serial import CGMappingSerial
from .scripts.cgmapping_residue import CGMappingResidue
from .scripts.fibo import FiboSphere

# THIS WILL WORK IF I DO
# from SMDAnalysis import *
#__all__ = ['common', 'compute']


### TO USE CYTHON DENSITY:
#cat > make.sh << EOF
##!/bin/bash
#cd compute
#python density_help_setup.py build_ext --inplace
#rm -rf build *.c *.html
#cd ..
#EOF
#
#bash make.sh
