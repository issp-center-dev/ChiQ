import numpy as np
from itertools import product
from collections import Counter
# from misc.point_group import PointGroup
from chiq.point_group_data import PointGroupData


def _check_num_irreps_symmetries(_irreps, _symmetries):
    # check duplicate
    def is_unique(array):
        return len(array) == len(set(array))

    assert is_unique([irp.name for irp in _irreps]), "Duplicate in names in irreps"
    assert is_unique([sym.name for sym in _symmetries]), "Duplicate in names in symmetries"

    # number of different kinds of symmetries (conjugacy class) = number of irreps
    assert len(_symmetries) == len(_irreps)

def _check_irreps(_irreps, _group_order):
    # sum dim(irp)^2 = group_order1
    print(" dims of irreps")
    sum_irp_dim_sqr = sum([irp.dim ** 2 for irp in _irreps])
    assert sum_irp_dim_sqr == _group_order, "Number of irreps inconsistent with the group order"

def _check_symmetries(_symmetries, _group_order):
    # sum sym.multiplicity = group_order1
    print(" multiplicities of symmetry operations")
    multiplicities = np.array([sym.multiplicity for sym in _symmetries])
    sum_multiplicities = np.sum(multiplicities)
    assert sum_multiplicities == _group_order, "Number of symmetries inconsistent with the group order"

def _check_character_table_size(_character_table, _irreps, _symmetries):
    # check size of character table
    print(" size of character table")
    assert len(_character_table) == len(_irreps), "size of the column"
    for chars in _character_table:
        assert len(chars) == len(_symmetries), "size of the row"

def _check_character_table_orthogonality(_character_table, _group_order, _func_product):
    # orthogonality
    print(" orthogonality of character table")
    for (i, chars1), (j, chars2) in product(enumerate(_character_table), repeat=2):
        prod = _func_product(chars1, chars2)
        if i == j:
            assert np.allclose(prod, _group_order), "(%d-th row) not normalized: %f" % (i, prod)
        else:
            assert np.allclose(prod, 0), "(%d, %d) rows not orthogonalized: %f" % (i, j, prod)

def _check_num_operators(_operators, _symmetries):
    # number of operators
    print(" number of operators")
    num_ops = {sym.name: sym.multiplicity for sym in _symmetries}
    print("  ", num_ops)
    num_ops_defined = Counter([op.name for op in _operators])
    print("  ", dict(num_ops_defined))
    assert num_ops_defined == num_ops, "# of operators defined is wrong"


def test_single_reps():
    for group in list(PointGroupData.keys()):
        print("\n=========== group=%s single reps =============" % repr(group))

        pg = PointGroupData[group](double_group=False)
        pg.use_single_reps()

        # pg.* = pg.*1
        _check_num_irreps_symmetries(pg.irreps, pg.symmetries)
        _check_irreps(pg.irreps, pg.group_order)
        _check_symmetries(pg.symmetries, pg.group_order)
        _check_character_table_size(pg.character_table, pg.irreps, pg.symmetries)
        _check_character_table_orthogonality(pg.character_table, pg.group_order, pg.inner_product)
        _check_num_operators(pg.operators, pg.symmetries)

def test_double_reps():
    for group in list(PointGroupData.keys()):
        print("\n=========== group=%s double reps =============" % repr(group))

        pg = PointGroupData[group](double_group=True)
        pg.use_double_reps()

        # pg.* = pg.*2           for [group_order, irreps, character_table]
        # pg.* = pg.*1 + pg.*2   for [symmetries, operators]

        # compare only double-group quantities
        _check_num_irreps_symmetries(pg.irreps2, pg.symmetries2)

        # compare single- + double-group quantities
        _check_irreps(pg.irreps1+pg.irreps2, pg.group_order2)

        # This check does not apply to double group
        # _check_symmetries(symmetries, group_order)

        _check_character_table_size(pg.character_table, pg.irreps, pg.symmetries)
        _check_character_table_orthogonality(pg.character_table, pg.group_order, pg.inner_product)
        _check_num_operators(pg.operators, pg.symmetries)

        # compare only double-group quantities
        _check_num_operators(pg.operators2, pg.symmetries2)
