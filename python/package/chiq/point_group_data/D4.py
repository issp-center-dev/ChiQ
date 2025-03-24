from .base import PointGroupDataBase

import math


class PointGroupData(PointGroupDataBase):

    def _define_pointgroup_data(self):

        self.group_order1 = 8

        IrpType = PointGroupDataBase.IrpType  # namedtuple('IrpType', ('name', 'dim'))
        self.irreps1 = [
            IrpType('A1', 1),
            IrpType('A2', 1),
            IrpType('B1', 1),
            IrpType('B2', 1),
            IrpType('E', 2),
        ]

        SymType = PointGroupDataBase.SymType  # namedtuple('SymType', ('name', 'multiplicity'))
        self.symmetries1 = [
            SymType('E', 1),
            SymType('C4', 2),
            SymType('C4^2', 1),
            SymType('C2a', 2),
            SymType('C2b', 2),
        ]

        self.character_table1 = [
            [1,  1,  1,  1,  1],
            [1,  1,  1, -1, -1],
            [1, -1,  1,  1, -1],
            [1, -1,  1, -1,  1],
            [2,  0, -2,  0,  0],
        ]

        OpType = PointGroupDataBase.OpType
        self.operators1 = [
            # E
            OpType('E', 'identity'),
            # 2 C4
            OpType('C4', 'rot', axis='z', angle='1/4'),
            OpType('C4', 'rot', axis='z', angle='-1/4'),
            # C4^2
            OpType('C4^2', 'rot', axis='z', angle='1/2'),
            # 2 C2'
            OpType('C2a', 'rot', axis='x', angle='1/2'),
            OpType('C2a', 'rot', axis='y', angle='1/2'),
            # 2 C2'
            OpType('C2b', 'rot', axis='x y', angle='1/2'),
            OpType('C2b', 'rot', axis='-x y', angle='1/2'),
        ]

        # for double group
        self.irreps2 =[
            IrpType('E1/2', 2),
            IrpType('E3/2', 2),
        ]

        self.character_table2 = [
            [2, math.sqrt(2), 0, 0, 0],
            [2, -math.sqrt(2), 0, 0, 0],
        ]
