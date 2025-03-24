

import numpy
from itertools import product
from collections import namedtuple, OrderedDict, Counter
import copy


class _OpType(object):
    def __init__(self, name, optype, axis=None, angle=None, inversion=False, double=False):
        self.name = name
        assert optype in ('identity', 'rot', 'mirror')
        self.optype = optype
        self.axis = axis
        self.angle = angle
        self.inversion = inversion
        self.double = double

        # mat should be set elsewhere
        self.mat = None
        self._mat_backup = None

    def backup(self):
        self._mat_backup = self.mat.copy()

    def restore(self):
        self.mat = self._mat_backup.copy()


class PointGroupDataBase(object):

    IrpType = namedtuple('IrpType', ('name', 'dim'))
    SymType = namedtuple('SymType', ('name', 'multiplicity'))
    # OpType = namedtuple('OpType', ('name', 'axis', 'angle', 'inversion'))
    OpType = _OpType

    def __init__(self, double_group=False):

        self._define_pointgroup_data()

        # check if necessary data have been defined
        for varname in ['group_order1', 'irreps1', 'symmetries1', 'operators1', 'character_table1']:
            assert hasattr(self, varname), "%s is not defined" % varname

        if double_group:
            # check if necessary data have been defined
            for varname in ['irreps2', 'character_table2']:
                assert hasattr(self, varname), "%s is not defined" % varname

            self._make_double_group()
        else:
            self.irreps2 = []
            self.character_table2 = []
            self.symmetries2 = []
            self.operators2 = []

        self.use_single_reps()

    def _define_pointgroup_data(self):
        """should be implemented in a subclass"""
        self.group_order1 = None
        self.irreps1 = None
        self.symmetries1 = None
        self.operators1 = None
        self.character_table1 = None
        raise NotImplementedError

    def use_single_reps(self):
        self.group_order = self.group_order1
        self.irreps = self.irreps1
        self.character_table = self.character_table1
        self.symmetries = self.symmetries1
        self.operators = self.operators1
        self._make_characters_dict()

    def use_double_reps(self):
        self.group_order = self.group_order2
        self.irreps = self.irreps2
        self.character_table = self.character_table2
        self.symmetries = self.symmetries1 + self.symmetries2
        self.operators = self.operators1 + self.operators2
        self._make_characters_dict()

    def inner_product(self, chars1, chars2):
        multiplicities = numpy.array([sym.multiplicity for sym in self.symmetries])
        assert len(chars1) == len(chars2) == len(multiplicities)
        return numpy.sum(numpy.array(chars1).conj() * numpy.array(chars2) * multiplicities)

    def make_inversion(self):
        self.group_order1 *= 2

        def make_inversion_irreps(_irreps):
            irreps_inv = []
            for irrep in _irreps:
                irreps_inv.append(PointGroupDataBase.IrpType(irrep.name + 'g', irrep.dim))
            for irrep in _irreps:
                irreps_inv.append(PointGroupDataBase.IrpType(irrep.name + 'u', irrep.dim))
            return irreps_inv

        self.irreps1 = make_inversion_irreps(self.irreps1)
        self.irreps2 = make_inversion_irreps(self.irreps2)

        symmetries_copy = copy.deepcopy(self.symmetries1)
        for sym in symmetries_copy:
            self.symmetries1.append(PointGroupDataBase.SymType('I' + sym.name, sym.multiplicity))

        def make_inversion_table(_table):
            # character_table_inv = copy.deepcopy(_character_table)
            table_inv = []
            for chars in _table:
                table_inv.append(chars + chars)
                # chars.extend(chars)
            for chars in _table:
                table_inv.append(chars + [-x for x in chars])
            return table_inv

        self.character_table1 = make_inversion_table(self.character_table1)
        self.character_table2 = make_inversion_table(self.character_table2)

        ops_inv = copy.deepcopy(self.operators1)
        for op in ops_inv:
            op.name = 'I' + op.name
            op.inversion = True
        self.operators1.extend(ops_inv)

    def _make_double_group(self):

        self.group_order2 = self.group_order1 * 2

        # find symmetries which has nonzero characters for double irreps
        flags_nonzero_char = [False] * len(self.symmetries1)
        for chars in self.character_table2:
            assert len(chars) == len(self.symmetries1)
            for i, c in enumerate(chars):
                if c != 0:
                    flags_nonzero_char[i] = True

        # list of symmetries to be extended
        syms_double = [sym for sym, flag_nonzero in zip(self.symmetries1, flags_nonzero_char) if flag_nonzero]
        syms_double_names = [sym.name for sym in syms_double]

        # symmetries for double reps
        self.symmetries2 = []
        for sym in syms_double:
            self.symmetries2.append(PointGroupDataBase.SymType('-' + sym.name, sym.multiplicity))

        # extend character table
        for chars in self.character_table2:
            chars_add = []
            for flag_nonzero, c in zip(flags_nonzero_char, chars):
                if flag_nonzero:
                    chars_add.append(-c)
            chars.extend(chars_add)

        # make symmetry operators for double reps
        self.operators2 = []
        for op in self.operators1:
            # print(op.name)
            if op.name in syms_double_names:
                op2 = copy.deepcopy(op)
                op2.name = '-' + op.name
                op2.double = True
                self.operators2.append(op2)

    def _make_characters_dict(self):
        # characters accessible by chars_dict[(irp_name, sym_name)]
        self.characters_dict = OrderedDict()
        sym_names = [sym.name for sym in self.symmetries]
        irp_names = [irp.name for irp in self.irreps]
        for (i, irp), (j, sym) in product(enumerate(irp_names), enumerate(sym_names)):
            self.characters_dict[(irp, sym)] = self.character_table[i][j]
        # print(self.chars_dict)
