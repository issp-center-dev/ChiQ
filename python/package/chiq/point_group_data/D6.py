from .base import PointGroupDataBase

import math


class PointGroupData(PointGroupDataBase):

    def _define_pointgroup_data(self):

        self.group_order1 = 12

        IrpType = PointGroupDataBase.IrpType  # namedtuple('IrpType', ('name', 'dim'))
        self.irreps1 = [
            IrpType('A1', 1),
            IrpType('A2', 1),
            IrpType('B1', 1),
            IrpType('B2', 1),
            IrpType('E1', 2),
            IrpType('E2', 2),
        ]

        SymType = PointGroupDataBase.SymType  # namedtuple('SymType', ('name', 'multiplicity'))
        self.symmetries1 = [
            SymType('E', 1),
            SymType('C6', 2),
            SymType('C3', 2),
            SymType('C2', 1),
            SymType('C2y', 3),
            SymType('C2x', 3),
        ]

        self.character_table1 = [
            [1,  1,  1,  1,  1,  1],
            [1,  1,  1,  1, -1, -1],
            [1, -1,  1, -1,  1, -1],
            [1, -1,  1, -1, -1,  1],
            [2,  1, -1, -2,  0,  0],
            [2, -1, -1,  2,  0,  0],
        ]

        OpType = PointGroupDataBase.OpType
        self.operators1 = [
            # E
            OpType('E', 'identity'),
            # 2 C6
            OpType('C6', 'rot', axis='z', angle='1/6'),
            OpType('C6', 'rot', axis='z', angle='-1/6'),
            # 2 C3
            OpType('C3', 'rot', axis='z', angle='1/3'),
            OpType('C3', 'rot', axis='z', angle='-1/3'),
            # C2
            OpType('C2', 'rot', axis='z', angle='1/2'),
            # 3 C2y
            OpType('C2y', 'rot', axis='y', angle='1/2'),
            OpType('C2y', 'rot', axis='1/12 xy', angle='1/2'),
            OpType('C2y', 'rot', axis='-1/12 xy', angle='1/2'),
            # 3 C2x
            OpType('C2x', 'rot', axis='x', angle='1/2'),
            OpType('C2x', 'rot', axis='1/6 xy', angle='1/2'),
            OpType('C2x', 'rot', axis='-1/6 xy', angle='1/2'),
        ]

        # for double group
        self.irreps2 =[
            IrpType('E1/2', 2),
            IrpType('E5/2', 2),
            IrpType('E3/2', 2),
        ]

        self.character_table2 = [
            [2, math.sqrt(3), 1, 0, 0, 0],
            [2, -math.sqrt(3), 1, 0, 0, 0],
            [2, 0, -2, 0, 0, 0],
        ]
