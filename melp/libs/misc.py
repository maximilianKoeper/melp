# ---------------------------------------------------------------------
#  miscellaneous functions
# ---------------------------------------------------------------------

import array


# ---------------------------------------------------------------------
#  searching index in arrays (fast)
# ---------------------------------------------------------------------
def index_finder(lst, item):
    """A generator function, if you might not need all the indices"""
    start = 0
    while True:
        try:
            start = lst.index(item, start)
            yield start
            start += 1
        except ValueError:
            break


def index_find_all(lst, item, results=None):
    """ If you want all the indices.
    Pass results=[] if you explicitly need a list,
    or anything that can .append(..)
    """
    if results is None:
        length = len(lst)
        results = (array.array('B') if length <= 2 ** 8 else
                   array.array('H') if length <= 2 ** 16 else
                   array.array('L') if length <= 2 ** 32 else
                   array.array('Q'))
    start = 0
    while True:
        try:
            start = lst.index(item, start)
            results.append(start)
            start += 1
        except ValueError:
            return results
