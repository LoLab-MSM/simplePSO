"""
Simple interface for particle swarm optimization

"""
from simplepso.pso import PSO

_MAJOR = 2
_MINOR = 2
_MICRO = 1
__version__ = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
__release__ = '%d.%d' % (_MAJOR, _MINOR)
__all__ = ['PSO']
