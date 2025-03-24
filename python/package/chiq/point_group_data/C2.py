from .base import PointGroupDataBase

import math


class PointGroupData(PointGroupDataBase):

    def _define_pointgroup_data(self):

        self.group_order1 = 2

        IrpType = PointGroupDataBase.IrpType  # namedtuple('IrpType', ('name', 'dim'))
        self.irreps1 = [
            IrpType('A', 1),
            IrpType('B', 1),
        ]

        SymType = PointGroupDataBase.SymType  # namedtuple('SymType', ('name', 'multiplicity'))
        self.symmetries1 = [
            SymType('E', 1),
            SymType('C2', 1),
        ]

        self.character_table1 = [
            [1,  1],
            [1, -1],
        ]

        OpType = PointGroupDataBase.OpType
        self.operators1 = [
            # E
            OpType('E', 'identity'),
            # C4^2
            OpType('C2', 'rot', axis='z', angle='1/2'),
        ]

        # for double group
        # self.irreps2 =[
        #     IrpType('E1/2', 2),
        # ]

        # self.character_table2 = [
        #     [2, 0],
        # ]

        self.irreps2 =[
            IrpType('E1/2(a)', 1),
            IrpType('E1/2(b)', 1),
        ]

        self.character_table2 = [
            [1, -1j],
            [1, 1j],
        ]
