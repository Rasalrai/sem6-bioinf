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

    return True  # differences not found


class FirstAlgorithm:
    def __init__(self, g: Graph):
        self.graph: Graph = g
        self.file: XmlFile = self.graph.xml

        # self.seq1_even = None
        # self.seq1_odd = None
        # self.seq2_even = None
        # self.seq2_odd = None

        # keeping track of our choices
        self.history1_even = []
        self.history1_odd = []
        self.history2_even = []
        self.history2_odd = []

    def execute(self):
        self.get_first_oligonucleotides()

        even = self.next_step(0)
        odd = self.next_step(1)

        while not (even or odd):
            even = self.next_step(0)
            odd = self.next_step(1)

        self.print_result()

    def get_first_oligonucleotides(self):
        first = "".join([self.file.start[i] if not i % 2 else "X" for i in range(len(self.file.start))])

        # there should be exactly and only one
        first_id = [a for a in range(len(self.graph.nodes1)) if self.graph.nodes1[a] == first][0]
        self.history1_even.append({"prev": -1, "choice": False, "chosen": first_id})

        second_start = "".join([self.file.start[i] if i % 2 else "X" for i in range(1, len(self.file.start))])
        sec_options = [on_id for on_id in range(len(self.graph.nodes1)) if
                       self.graph.nodes1[on_id].startswith(second_start)]

        second_id = sec_options[0]
        if len(sec_options) > 1:
            # save it as the first choice
            self.history1_odd.append({"prev": -1, "chosen": second_id, "choice": True, "options": sec_options[1:]})
        else:
            self.history1_odd.append({"prev": -1, "chosen": second_id, "choice": False})

        print(f"sequence 1 starting with: {first_id} and {second_id}\n{self.file.start}")
        print(self.graph.nodes1[first_id], self.graph.nodes1[second_id], sep="\n·")

        # print(self.history1_even, self.history1_odd, sep="\n")
        # print(f"\nnext options (s1):\n{self.get_next_nodes(first_id, 1)}\t{self.get_next_nodes(second_id, 1)}\n\n")

    def get_next_nodes(self, node_no, spectr_no):
        """get list of direct successors of node `node_no` in graph im1 or im2"""
        if spectr_no == 1:
            return [p for p in range(len(self.graph.im1[node_no])) if self.graph.im1[node_no][p]]
        elif spectr_no == 2:
            return [p for p in range(len(self.graph.im2[node_no])) if self.graph.im2[node_no][p]]

    def next_step(self, parity: int):
        """
        append the next node to the list
        :param parity: 0 for even, 1 for odd
        """

        history = self.history1_odd if parity % 2 else self.history1_even
        entry = {"prev": history[-1]["chosen"]}

        options = self.get_next_nodes(entry["prev"], 1)

        if len(options) == 0:
            print(f"-- cannot continue; no nodes to choose from ({'odd' if parity % 2 else 'even'}) --")
            return -1

        entry["chosen"] = options[0]

        if len(options) > 1:
            entry["choice"] = True
            entry["options"] = options[1:]
        else:
            entry["choice"] = False

        history.append(entry)
        return 0

    def print_result(self):
        # TODO it should work also if even and odd lengths are different

        print("\n\n=== RESULT ===\n")
        spaces = 0
        sequence = self.file.start[:-1]

        for e, o in zip(self.history1_even, self.history1_odd):
            print(" "*spaces, self.graph.nodes1[e["chosen"]], sep="")
            spaces += 1
            print(" "*spaces, self.graph.nodes1[o["chosen"]], sep="")
            spaces += 1

            sequence += self.graph.nodes1[e["chosen"]][-1]
            sequence += self.graph.nodes1[o["chosen"]][-1]

        print(f"\n{sequence}\n  length: {spaces + len(self.graph.nodes1[0]) - 1}")
