from collections import namedtuple
from rdkit import Chem as Chem
from rdkit.Chem import Crippen as Crippen, MolSurf as MolSurf
from rdkit.Chem.ChemUtils.DescriptorUtilities import setDescriptorVersion as setDescriptorVersion
from typing import Any, Optional

QEDproperties = namedtuple('QEDproperties', 'MW,ALOGP,HBA,HBD,PSA,ROTB,AROM,ALERTS')

ADSparameter = namedtuple('ADSparameter', 'A,B,C,D,E,F,DMAX')
WEIGHT_MAX: Any
WEIGHT_MEAN: Any
WEIGHT_NONE: Any
AliphaticRings: Any
AcceptorSmarts: Any
Acceptors: Any
StructuralAlertSmarts: Any
StructuralAlerts: Any
adsParameters: Any

def ads(x: Any, adsParameter: Any): ...
def properties(mol: Any): ...
def qed(mol: Any, w: Any = ..., qedProperties: Optional[Any] = ...): ...
def weights_max(mol: Any): ...
def weights_mean(mol: Any): ...
def weights_none(mol: Any): ...
def default(mol: Any): ...
