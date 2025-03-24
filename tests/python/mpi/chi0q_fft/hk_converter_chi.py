from dcore.converters.hk import *
import numpy as np

class HkConverterChi(HkConverter):

    # __init__ is not changed from that in the base class

    # def convert_dft_input_chi(self, div, kpoint_ibz, k2i):
    #     subgrp = self.dft_subgrp + '_chi'
    #     ar = HDFArchive(self.hdf_file, 'a')
    #     if not (subgrp in ar):
    #         ar.create_group(subgrp)
    #     ar[subgrp]['div'] = np.array(div)
    #     ar[subgrp]['k_points_ibz'] = np.array(kpoint_ibz)
    #     ar[subgrp]['k_to_ik'] = k2i
    #     del ar

    def convert_dft_input_chi(self, div, chi_subgroup='dft_input_chi'):
        # subgrp = self.dft_subgrp + '_chi'
        subgrp = chi_subgroup
        ar = HDFArchive(self.hdf_file, 'a')
        if not (subgrp in ar):
            ar.create_group(subgrp)
        ar[subgrp]['div'] = np.array(div)
        del ar
