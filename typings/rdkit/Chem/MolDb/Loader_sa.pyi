from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem, Crippen as Crippen, Descriptors as Descriptors, Lipinski as Lipinski
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from sqlalchemy import DateTime as DateTime, ForeignKey as ForeignKey, MetaData as MetaData, Table as Table
from sqlalchemy.orm import backref as backref, mapper as mapper, relation as relation
from typing import Any, Optional

decBase: Any

class Compound(decBase):
    __tablename__: str = ...
    guid: Any = ...
    molpkl: Any = ...

def RegisterSchema(dbUrl: Any, echo: bool = ...): ...
ConnectToSchema = RegisterSchema
logger: Any

def ProcessMol(session: Any, mol: Any, globalProps: Any, nDone: Any, nameProp: str = ..., nameCol: str = ..., redraw: bool = ..., keepHs: bool = ..., skipProps: bool = ..., addComputedProps: bool = ..., skipSmiles: bool = ...): ...
def LoadDb(suppl: Any, dbName: Any, nameProp: str = ..., nameCol: str = ..., silent: bool = ..., redraw: bool = ..., errorsTo: Optional[Any] = ..., keepHs: bool = ..., defaultVal: str = ..., skipProps: bool = ..., regName: str = ..., skipSmiles: bool = ..., maxRowsCached: int = ..., uniqNames: bool = ..., addComputedProps: bool = ..., lazySupplier: bool = ..., numForPropScan: int = ..., startAnew: bool = ...) -> None: ...
