from .base import PointGroupDataBase

import math


class PointGroupData(PointGroupDataBase):

    def _define_pointgroup_data(self):

        self.group_order1 = 24

        IrpType = PointGroupDataBase.IrpType  # namedtuple('IrpType', ('name', 'dim'))
        self.irreps1 = [
            IrpType('A1', 1),
            IrpType('A2', 1),
            IrpType('E', 2),
            IrpType('T1', 3),
            IrpType('T2', 3),
        ]

        SymType = PointGroupDataBase.SymType  # namedtuple('SymType', ('name', 'multiplicity'))
        self.symmetries1 = [
            SymType('E', 1),
            SymType('C4', 6),
            SymType('C4^2', 3),
            SymType('C2a', 6),
            SymType('C3', 8),
        ]

        self.character_table1 = [
            [1,  1,  1,  1,  1],
            [1, -1,  1, -1,  1],
            [2,  0,  2,  0, -1],
            [3,  1, -1, -1,  0],
            [3, -1, -1,  1,  0],
        ]

        OpType = PointGroupDataBase.OpType
        self.operators1 = [
            # E
            OpType('E', 'identity'),
            # 6 C4
            OpType('C4', 'rot', axis='z', angle='1/4'),
            OpType('C4', 'rot', axis='z', angle='-1/4'),
            OpType('C4', 'rot', axis='x', angle='1/4'),
            OpType('C4', 'rot', axis='x', angle='-1/4'),
            OpType('C4', 'rot', axis='y', angle='1/4'),
            OpType('C4', 'rot', axis='y', angle='-1/4'),
            # 3 C4^2
            OpType('C4^2', 'rot', axis='z', angle='1/2'),
            OpType('C4^2', 'rot', axis='x', angle='1/2'),
            OpType('C4^2', 'rot', axis='y', angle='1/2'),
            # 6 C2a
            OpType('C2a', 'rot', axis='x y', angle='1/2'),
            OpType('C2a', 'rot', axis='-x y', angle='1/2'),
            OpType('C2a', 'rot', axis='y z', angle='1/2'),
            OpType('C2a', 'rot', axis='-y z', angle='1/2'),
            OpType('C2a', 'rot', axis='z x', angle='1/2'),
            OpType('C2a', 'rot', axis='-z x', angle='1/2'),
            # 8 C3
            OpType('C3', 'rot', axis='x y z', angle='1/3'),
            OpType('C3', 'rot', axis='x y z', angle='-1/3'),
            OpType('C3', 'rot', axis='-x y z', angle='1/3'),
            OpType('C3', 'rot', axis='-x y z', angle='-1/3'),
            OpType('C3', 'rot', axis='x -y z', angle='1/3'),
            OpType('C3', 'rot', axis='x -y z', angle='-1/3'),
            OpType('C3', 'rot', axis='x y -z', angle='1/3'),
            OpType('C3', 'rot', axis='x y -z', angle='-1/3'),
        ]

        # for double group
        self.irreps2 =[
            IrpType('E1/2', 2),
            IrpType('E5/2', 2),
            IrpType('G3/2', 4),
        ]

        self.character_table2 = [
            [2, math.sqrt(2), 0, 0, 1],
            [2, -math.sqrt(2), 0, 0, 1],
            [4, 0, 0, 0, -1],
        ]
