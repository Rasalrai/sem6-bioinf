import random

from Graph import Graph
from XmlFile import XmlFile


def match_sequences(first, second, shift=0):
    """
    e.g. with shift=2
    AXAXGXCXT   first
    ..AXGXCXTG  second

    only works for shift >= 0
    """
    f_id = 0
    s_id = shift
    while f_id < len(first) and s_id < len(second):
        # if the letters are different and they are not X, the sequences don't match
        if first[f_id] != second[s_id] and "X" not in (first[f_id], second[s_id]):
            return False
        f_id += 1
        s_id += 1

    return True     # differences not found


class FirstAlgorithm:
    def __init__(self, f: XmlFile, g: Graph):
        self.file = f
        self.graph = g
        self.seq1 = None    # odd
        self.seq2 = None    # even
        self.history1 = []
        self.history2 = []

    def get_first_oligonucleotides(self):
        first = [self.file.start[i] if not i % 2 else "X" for i in range(len(self.file.start))]

        second_start = "".join([self.file.start[i] if i % 2 else "X" for i in range(1, len(self.file.start))])
        sec_options = [on_id for on_id in range(len(self.file.s1)) if self.file.s1[on_id].startswith(second_start)]

        second = sec_options[0]
        if len(sec_options) > 1:
            # save it as the first choice
            self.history2.append({"prev": -1, "chosen": second, "options": sec_options})

        self.seq1 = "".join(first)
        self.seq2 = self.file.s1[second]

    def get_next_nodes(self, node_no, graph_no):
        """get list of direct successors of node `node_no` in graph im1 or im2"""
        if graph_no == 1:
            return [p for p in range(len(self.graph.im1[node_no])) if self.graph.im1[node_no][p]]
        elif graph_no == 2:
            return [p for p in range(len(self.graph.im2[node_no])) if self.graph.im2[node_no][p]]


