from .base import PointGroupDataBase

import math


class PointGroupData(PointGroupDataBase):

    def _define_pointgroup_data(self):

        self.group_order1 = 1

        IrpType = PointGroupDataBase.IrpType  # namedtuple('IrpType', ('name', 'dim'))
        self.irreps1 = [
            IrpType('A', 1),
        ]

        SymType = PointGroupDataBase.SymType  # namedtuple('SymType', ('name', 'multiplicity'))
        self.symmetries1 = [
            SymType('E', 1),
        ]

        self.character_table1 = [
            [1,],
        ]

        OpType = PointGroupDataBase.OpType
        self.operators1 = [
            # E
            OpType('E', 'identity'),
        ]

        # for double group
        self.irreps2 =[
            IrpType('B1/2', 1),
        ]

        self.character_table2 = [
            [1,],
        ]
