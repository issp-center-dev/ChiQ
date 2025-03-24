class Communicator:
    def __init__(self):
        pass

    def gather(self, data, root=0):
        return [data]

    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1

    def barrier(self):
        pass

    def bcast(self, data, root=0):
        return data

    def Abort(self, errorcode=0):
        import sys
        sys.exit(errorcode)

COMM_WORLD = Communicator()
