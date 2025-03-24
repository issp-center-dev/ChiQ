from collections import namedtuple

IrpType = namedtuple('IrpType', ('name', 'dim'))
SymType = namedtuple('SymType', ('name', 'multiplicity'))
OpType = namedtuple('OpType', ('name', 'axis', 'angle', 'inversion'))
