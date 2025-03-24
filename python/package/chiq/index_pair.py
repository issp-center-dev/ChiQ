from itertools import product


class IndexPair:
    """
    generator of product of indices

    Ex.)
    indices = ['up','down']
    -->
    namelist = ['up-up', 'up-down', 'down-up', 'down-down']
    pairlist = [('up', 'up'), ('up', 'down'), ('down', 'up'), ('down', 'down')]

    get_name(('up','down')) --> 'up-down'
    get_index(('up','down')) --> 1
    """

    def __init__(self, indices, separator="-", convert_to_int=False, only_diagonal=False):
        """
        indices : list<string> or list<int>
        convert_to_int: if True, indices are converted to zero-based integer
        only_diagonal: if True, only diagonal components are taken, e.g., ('up', 'up') and ('down', 'down')
        """
        self.indices = indices
        # self.size = len(indices)**2

        self.namelist = []
        self.pairlist = []
        self.__name_pair_list = []  # for iterator
        self.__converter_name = {}
        self.__converter_index = {}
        count = 0
        for (i1,s1), (i2,s2) in product(enumerate(indices), repeat=2):
            if only_diagonal and i1 != i2:
                continue
            name = str(s1) + separator + str(s2)
            pair = (i1,i2) if convert_to_int else (s1,s2)
            self.namelist.append(name)
            self.pairlist.append(pair)
            self.__name_pair_list.append((name,pair))
            self.__converter_name[pair] = name
            self.__converter_index[pair] = count
            count += 1

        assert count == len(self.namelist)
        self.size = count

    def __iter__(self):
        return iter(self.__name_pair_list)

    def get_name(self, index1, index2):
        return self.__converter_name[(index1, index2)]

    def get_index(self, index1, index2):
        return self.__converter_index[(index1, index2)]


class IndexPair2:
    """
    generator of product of indices

    Ex.)
    indices1 = [0,1]  # two-sublattice
    indices2 = ['up','down']
    -->
    namelist = ['0-up-0-up', '0-up-0-down', '0-down-0-up', '0-down-0-down',
                '0-up-1-up', '0-up-1-down', '0-down-1-up', '0-down-1-down',
                '1-up-0-up', '1-up-0-down', '1-down-0-up', '1-down-0-down',
                '1-up-1-up', '1-up-1-down', '1-down-1-up', '1-down-1-down',
    pairlist = [(0,'up',0,'up'), (0,'up',0,'down'), (0,'down',0,'up'), (0,'down',0,'down')
                (0,'up',1,'up'), (0,'up',1,'down'), (0,'down',1,'up'), (0,'down',1,'down')
                (1,'up',0,'up'), (1,'up',0,'down'), (1,'down',0,'up'), (1,'down',0,'down')
                (1,'up',1,'up'), (1,'up',1,'down'), (1,'down',1,'up'), (1,'down',1,'down')]

    get_name((0, 'up', 1, 'down')) --> '0-up-1-down'
    get_index((0, 'up', 1, 'down')) --> 5
    """

    def __init__(self, indices1, indices2, separator="-", convert_to_int1=False, convert_to_int2=False, only_diagonal1=False, only_diagonal2=False, indices_order='1212'):
        """
        indices : list<string> or list<int>
        convert_to_int: if True, indices are converted to zero-based integer
        only_diagonal: if True, only diagonal components are taken, e.g., ('up', 'up') and ('down', 'down')
        """
        self.indices1 = indices1
        self.indices2 = indices2
        # self.size = (len(indices1)*len(indices2))**2

        self.namelist = []
        self.pairlist = []
        self.__name_pair_list = []  # for iterator
        self.__converter_name = {}
        self.__converter_index = {}
        count = 0
        for (i1,s1), (i2,s2), (j1,t1), (j2,t2) in product(enumerate(indices1), enumerate(indices1), enumerate(indices2), enumerate(indices2)):
            if only_diagonal1 and i1 != i2:
                continue
            if only_diagonal2 and j1 != j2:
                continue
            name = str(s1) + separator + str(t1) + separator + str(s2) + separator + str(t2)
            # pair = (s1,t1,s2,t2)
            pair = (i1 if convert_to_int1 else s1,
                    j1 if convert_to_int2 else t1,
                    i2 if convert_to_int1 else s2,
                    j2 if convert_to_int2 else t2)
            self.namelist.append(name)
            self.pairlist.append(pair)
            self.__name_pair_list.append((name,pair))
            self.__converter_name[pair] = name
            self.__converter_index[pair] = count
            count += 1

        assert count == len(self.namelist)
        self.size = count

    def __iter__(self):
        return iter(self.__name_pair_list)

    def get_name(self, index1, index2, index3, index4):
        return self.__converter_name[(index1, index2, index3, index4)]

    def get_index(self, index1, index2, index3, index4):
        return self.__converter_index[(index1, index2, index3, index4)]
