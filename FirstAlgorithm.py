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
        self.on_len = len(self.graph.nodes1[0])     # oligonucleotide length

        # keeping track of the choices/options
        # probe 1
        self.history1_even = []
        self.history1_odd = []
        # probe 2
        self.history2_even = []
        self.history2_odd = []

    def execute(self):
        self.get_first_oligonucleotides()

        # even_length = ceil(total _len / 2)
        # odd_length = floor(total _len / 2)
        odd_len = (self.file.seq_len - self.on_len + 1) // 2
        evn_len = self.file.seq_len - self.on_len + 1 - odd_len

        even_chk = self.next_step(0)
        while len(self.history1_even) < evn_len:
            if even_chk:
                self.go_back_dfs(0)
            even_chk = self.next_step(0)

        odd_chk = self.next_step(1)
        while len(self.history1_odd) < odd_len:
            if odd_chk:
                self.go_back_dfs(1)
            odd_chk = self.next_step(1)

        # # while not (even or odd):
        # while not (even and odd):
        #     even = self.next_step(0)
        #     odd = self.next_step(1)

        self.print_result()

    def get_first_oligonucleotides(self):
        first = "".join([self.file.start[i] if not i % 2 else "X" for i in range(len(self.file.start))])

        # there should be exactly and only one
        first_id = [a for a in range(len(self.graph.nodes1)) if self.graph.nodes1[a] == first][0]
        self.history1_even.append({"prev": -1, "fork": False, "chosen": first_id})

        second_start = "".join([self.file.start[i] if i % 2 else "X" for i in range(1, len(self.file.start))])
        sec_options = [on_id for on_id in range(len(self.graph.nodes1)) if
                       self.graph.nodes1[on_id].startswith(second_start)]

        second_id = sec_options[0]
        if len(sec_options) > 1:
            # save it as the first choice
            self.history1_odd.append({"prev": -1, "chosen": second_id, "fork": True, "options": sec_options[1:]})
        else:
            self.history1_odd.append({"prev": -1, "chosen": second_id, "fork": False})

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
            print(f"- Cannot continue; no nodes to choose from ({'odd' if parity % 2 else 'even'}) --")
            return -1

        entry["chosen"] = options[0]

        if len(options) > 1:
            entry["fork"] = True
            entry["options"] = options[1:]
        else:
            entry["fork"] = False

        history.append(entry)
        return 0

    def go_back_dfs(self, parity):
        """ return to the last encountered fork """
        # todo go back in the probe 2 history too
        history = self.history1_odd if parity % 2 else self.history1_even
        # get the last spot with choices
        for i in reversed(range(len(history))):
            if history[i]["fork"]:
                # return to this point
                print(f"Return sequence {'odd' if parity % 2 else 'even'} from {len(history)} to {i}")
                history = history[:i+1]
                history[i]["chosen"] = history[i]["options"].pop(0)
                if len(history[i]["options"]) == 0:
                    history[i]["fork"] = False

                # prevent incorrect copies of the array (unexpected behaviour)
                if parity % 2:
                    self.history1_odd = history
                else:
                    self.history1_even = history
                return
        print("Did not find a fork to return to!")

        # return to this spot

    def print_result(self):
        print("\n\n=== RESULTS ===\n")
        spaces = 0
        sequence = self.file.start[:-1]

        for e, o in zip(self.history1_even, self.history1_odd):
            print(" "*spaces, self.graph.nodes1[e["chosen"]], sep="")
            spaces += 1
            print(" "*spaces, self.graph.nodes1[o["chosen"]], sep="")
            spaces += 1

            sequence += self.graph.nodes1[e["chosen"]][-1]
            sequence += self.graph.nodes1[o["chosen"]][-1]

        if len(self.history1_even) > len(self.history1_odd):
            e = self.history1_even[-1]
            print(" " * spaces, self.graph.nodes1[e["chosen"]], sep="")
            spaces += 1
            sequence += self.graph.nodes1[e["chosen"]][-1]

        print(f"\n{sequence}\n  length: {spaces + self.on_len - 1}")

    def print_single(self, parity):
        """ for debugging purposes - it's not called anywhere in the algorithm """
        history = self.history1_odd if parity % 2 else self.history1_even

        print("\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        spaces = 0
        for item in history:
            print("  " * spaces, self.graph.nodes1[item["chosen"]], " -F-" if item["fork"] else "", sep="")
            spaces += 1
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n\n")
