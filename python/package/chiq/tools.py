import time


class WallTime(object):
    def __init__(self, str_walltime="\nWallTime:", str_time="\nTime:"):
        self._t = [time.time()] * 2
        self._str_walltime = str_walltime
        self._str_time = str_time

    def print_wall_time(self, comment):
        t_now = time.time()
        print("{} {:.1f} sec ({:.1f} sec for {})".format(self._str_walltime, t_now - self._t[0], t_now - self._t[1], comment))
        self._t[1] = t_now

    def print_time(self, comment):
        t_now = time.time()
        print("{} {:.1f} sec for {}".format(self._str_time, t_now - self._t[1], comment))
        self._t[1] = t_now
