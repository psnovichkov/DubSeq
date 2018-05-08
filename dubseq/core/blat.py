class BlatRecord:

    __type_converter = {
        's': lambda value: value,
        'd': lambda value: int(value),
        'a': lambda value: [int(x) for x in value.split(',')[:-1]]
    }

    def __init__(self):
        self.__match = None
        self.__mismatch = None
        self.__repmatch = None
        self.__ns = None
        self.__q_gapcount = None
        self.__q_gapbases = None
        self.__t_gapcount = None
        self.__t_gapbases = None
        self.__strand = None
        self.__q_name = None
        self.__q_size = None
        self.__q_start = None
        self.__q_end = None
        self.__t_name = None
        self.__t_size = None
        self.__t_start = None
        self.__t_end = None
        self.__block_count = None
        self.__block_sizes = None
        self.__q_starts = None
        self.__t_starts = None

    def __call__(self, line):
        vals = line.split('\t')
        (
            self.__match, self.__mismatch, self.__repmatch, self.__ns,
            self.__q_gapcount, self.__q_gapbases,
            self.__t_gapcount, self.__t_gapbases,
            self.__strand,
            self.__q_name, self.__q_size, self.__q_start, self.__q_end,
            self.__t_name, self.__t_size, self.__t_start, self.__t_end,
            self.__block_count, self.__block_sizes, self.__q_starts, self.__t_starts
        ) = (
            BlatRecord.__type_converter[x](vals[i])
            for i, x in enumerate('ddddddddssdddsddddaaa')
        )

    @property
    def match(self): return self.__match

    @property
    def mismatch(self): return self.__mismatch

    @property
    def repmatch(self): return self.__repmatch

    @property
    def ns(self): return self.__ns

    @property
    def q_gapcount(self): return self.__q_gapcount

    @property
    def q_gapbases(self): return self.__q_gapbases

    @property
    def t_gapcount(self): return self.__t_gapcount

    @property
    def t_gapbases(self): return self.__t_gapbases

    @property
    def strand(self): return self.__strand

    @property
    def q_name(self): return self.__q_name

    @property
    def q_size(self): return self.__q_size

    @property
    def q_start(self): return self.__q_start

    @property
    def q_end(self): return self.__q_end

    @property
    def t_name(self): return self.__t_name

    @property
    def t_size(self): return self.__t_size

    @property
    def t_start(self): return self.__t_start

    @property
    def t_end(self): return self.__t_end

    @property
    def block_count(self): return self.__block_count

    @property
    def block_sizes(self): return self.__block_sizes

    @property
    def q_starts(self): return self.__q_starts

    @property
    def t_starts(self): return self.__t_starts


class BlatReader:

    def __init__(self, file_name):
        self.__file_name = file_name
        self.__file = open(self.__file_name, 'r')

        # position a file cursor to the first record
        while True:
            line = self.__file.readline()
            if line.startswith('----'):
                break

    @property
    def file_name(self): return self.__file_name

    def next_record(self, blat_record):
        line = self.__file.readline()
        if not line:
            return False

        blat_record(line)
        return True

    def close(self):
        self.__file.close()
