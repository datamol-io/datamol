from rdkit.RDPaths import *
from pyPgSQL import PgSQL as PgSQL
from pysqlite2 import dbapi2 as dbapi2
from typing import Any

RDBaseDir: Any
RDCodeDir: Any
RDDataDir: Any
RDDocsDir: Any
RDDemoDir: Any
RDBinDir: Any
RDProjDir: Any
RDContribDir: Any
rpcTestPort: int
pythonTestCommand: str
defaultDBUser: str
defaultDBPassword: str

class ObsoleteCodeError(Exception): ...
class UnimplementedCodeError(Exception): ...

pythonExe: Any
usePgSQL: bool
useSqlLite: bool
RDTestDatabase: str
RDDataDatabase: str
molViewer: Any
