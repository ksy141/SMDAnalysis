import os
from pathlib import Path

_PKG = Path(__file__).parent

for _name, _var in [
    ('FF/CHARMM/toppar', 'CHARMM_TOPPAR'),
    ('mapping',          'SMDA_MAPPING'),
    ('profiles',         'SMDA_PROFILES'),
]:
    _path = _PKG / _name
    if _path.is_dir():
        os.environ.setdefault(_var, str(_path))
