from .share import *

group_order = 24

inversion = False

# IrpType = namedtuple('IrpType', ('name', 'dim'))
irreps = []
irreps.append(IrpType('A1', 1))
irreps.append(IrpType('A2', 1))
irreps.append(IrpType('E', 2))
irreps.append(IrpType('T1', 3))
irreps.append(IrpType('T2', 3))

# SymType = namedtuple('SymType', ('name', 'multiplicity'))
symmetries = []
symmetries.append(SymType('E', 1))
symmetries.append(SymType('C3', 8))
symmetries.append(SymType('C2', 6))
symmetries.append(SymType('C4', 6))
symmetries.append(SymType('C4^2', 3))

character_table = (
    (1,  1,  1,  1,  1),
    (1,  1, -1, -1,  1),
    (2, -1,  0,  0,  2),
    (3,  0, -1,  1, -1),
    (3,  0,  1, -1, -1),
)

# OpType = namedtuple('OpType', ('name', 'axis', 'angle', 'inversion'))
operators = []
# E
operators.append(OpType('E', '', '0', False))
# 8 C3
operators.append(OpType('C3', 'x y z', '1/3', False))
operators.append(OpType('C3', 'x y z', '2/3', False))
operators.append(OpType('C3', '-x y z', '1/3', False))
operators.append(OpType('C3', '-x y z', '2/3', False))
operators.append(OpType('C3', 'x -y z', '1/3', False))
operators.append(OpType('C3', 'x -y z', '2/3', False))
operators.append(OpType('C3', 'x y -z', '1/3', False))
operators.append(OpType('C3', 'x y -z', '2/3', False))
# 6 C2
operators.append(OpType('C2', 'x y', '1/2', False))
operators.append(OpType('C2', '-x y', '1/2', False))
operators.append(OpType('C2', 'y z', '1/2', False))
operators.append(OpType('C2', '-y z', '1/2', False))
operators.append(OpType('C2', 'z x', '1/2', False))
operators.append(OpType('C2', '-z x', '1/2', False))
# 6 C4
operators.append(OpType('C4', 'x', '1/4', False))
operators.append(OpType('C4', 'x', '3/4', False))
operators.append(OpType('C4', 'y', '1/4', False))
operators.append(OpType('C4', 'y', '3/4', False))
operators.append(OpType('C4', 'z', '1/4', False))
operators.append(OpType('C4', 'z', '3/4', False))
# 3 C4^2
operators.append(OpType('C4^2', 'x', '1/2', False))
operators.append(OpType('C4^2', 'y', '1/2', False))
operators.append(OpType('C4^2', 'z', '1/2', False))
