from .O import *


class PointGroupData(PointGroupData):

    def _define_pointgroup_data(self):
        super(PointGroupData, self)._define_pointgroup_data()
        self.make_inversion()
