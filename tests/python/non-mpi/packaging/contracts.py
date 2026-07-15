CANONICAL_COMMANDS = (
    "chiq_main", "chiq_post", "chiq_fft", "gen_qpath", "gen_allq",
    "calc_Iq", "calc_Iq_scl", "plot_chiq_path", "plot_Ir",
    "eigenvec_viewer", "dcore_chiq",
)

DEPRECATED_COMMANDS = tuple(name + ".py" for name in CANONICAL_COMMANDS)

SHIM_MODULES = (
    "bse_toml", "g2scl_core", "h5bse", "index_pair", "matrix_dict",
    "mpi", "point_group", "sumk_dft_chi", "tools",
)

POINT_GROUP_MODULES = ("C1", "C2", "D3", "D4", "D6", "O", "Oh", "base")

NATIVE_SUFFIXES = (".so", ".pyd", ".dylib")
