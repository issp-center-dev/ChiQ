from chiq import __version__

VERSION_MESSAGE = "ChiQ version %s" % __version__


def add_version_argument(parser):
    parser.add_argument("--version", action="version", version=VERSION_MESSAGE)
