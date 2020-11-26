from pyPgSQL.PgSQL import *
from rdkit import RDConfig as RDConfig
from typing import Any

sqlTextTypes: Any
sqlIntTypes: Any
sqlFloatTypes: Any
sqlBinTypes: Any
getTablesSql: str
getTablesAndViewsSql: str
getDbSql: str
fileWildcard: Any
placeHolder: str
binaryTypeName: str
binaryHolder = PgBytea
RDTestDatabase: str
dbFileWildcard: str
fileWildcard = dbFileWildcard
binaryHolder = memoryview

def connect(x: Any, *args: Any): ...
