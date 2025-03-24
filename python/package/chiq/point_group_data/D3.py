from .base import PointGroupDataBase

import math


class PointGroupData(PointGroupDataBase):

    def _define_pointgroup_data(self):

        self.group_order1 = 6

        IrpType = PointGroupDataBase.IrpType  # namedtuple('IrpType', ('name', 'dim'))
        self.irreps1 = [
            IrpType('A1', 1),
            IrpType('A2', 1),
            IrpType('E', 2),
        ]

        SymType = PointGroupDataBase.SymType  # namedtuple('SymType', ('name', 'multiplicity'))
        self.symmetries1 = [
            SymType('E', 1),
            SymType('C3', 2),
            SymType('C2x', 3),
        ]

        self.character_table1 = [
            [1,  1,  1],
            [1,  1, -1],
            [2, -1,  0],
        ]

        OpType = PointGroupDataBase.OpType
        self.operators1 = [
            # E
            OpType('E', 'identity'),
            # 2 C3
            OpType('C3', 'rot', axis='z', angle='1/3'),
            OpType('C3', 'rot', axis='z', angle='-1/3'),
            # 3 C2x
            OpType('C2x', 'rot', axis='x', angle='1/2'),
            OpType('C2x', 'rot', axis='1/6 xy', angle='1/2'),
            OpType('C2x', 'rot', axis='-1/6 xy', angle='1/2'),
        ]

        # for double group
        # self.irreps2 =[
        #     IrpType('E1/2', 2),
        #     IrpType('E3/2', 2),
        # ]

        # self.character_table2 = [
        #     [2, 1, 0],
        #     [2, -2, 0],
        # ]

        self.irreps2 =[
            IrpType('E1/2', 2),
            IrpType('E3/2(a)', 1),
            IrpType('E3/2(b)', 1),
        ]

        self.character_table2 = [
            [2, 1, 0],
            [1, -1, 1j],
            [1, -1, -1j],
        ]
