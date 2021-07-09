import copy
import random

from Graph import Graph
from XmlFile import XmlFile


class HeuristicAlgorithm:
    def __init__(self, g: Graph):
        self.graph: Graph = g
        self.file: XmlFile = self.graph.xml
        self.on1_len = len(self.graph.nodes1[0])  # oligonucleotide length

        # keeping track of the choices/options
        # probe 1
        self.history1_even = []
        self.history1_odd = []
        # probe 2
        self.history2_even = []
        self.history2_odd = []
        self.final_sequence = ""

    def build_paths(self):
        odd_len = (self.file.seq_len - self.on1_len + 1) // 2
        evn_len = self.file.seq_len - self.on1_len + 1 - odd_len

        seq2_evn_chk = 0
        seq2_odd_chk = 0

        evn_chk = 0
        odd_chk = 0

        # if s2 check returns !=0, then don't check until the next fork return
        # while len(self.history1_even) + len(self.history1_odd) < self.file.seq_len - self.on1_len + 1:
        while len(self.history1_even) <= len(self.history1_odd) and len(self.history1_even) < evn_len:
            if evn_chk:
                self.go_back_dfs(0, True)
            evn_chk = self.next_step(0)

        while len(self.history1_odd) < len(self.history1_even) and len(self.history1_odd) < odd_len:
            if odd_chk:
                self.go_back_dfs(1, True)
            odd_chk = self.next_step(1)

        # check if it's there was no error since last fork and if there's enough data from seq1
        # to check with seq2
        while not seq2_evn_chk and len(self.history2_even) <= len(self.history1_even):
            seq2_evn_chk = self.next_s2(0)

        if seq2_evn_chk:
            return 2

        while not seq2_odd_chk and len(self.history2_odd) < len(self.history1_odd):
            seq2_odd_chk = self.next_s2(1)

        if seq2_odd_chk:
            return 1

        return 0

    def execute(self):
        self.get_first_oligonucleotides()

        # even_length = ceil(total _len / 2)
        # odd_length = floor(total _len / 2)
        total_len = self.file.seq_len - self.on1_len + 1
        odd_len = total_len // 2
        evn_len = total_len - odd_len

        # if s2 check returns !=0, then don't check until the next fork return
        while len(self.history1_even) + len(self.history1_odd) < total_len or\
                len(self.history2_even) + len(self.history2_odd) <= total_len:
            if self.build_paths():
                rand_seq = random.randint(0, 1)
                # change a random path, with returns possible
                # if return not possible, return the other sequence
                if self.go_back_dfs(rand_seq, True):
                    self.go_back_dfs(rand_seq+1, True)

        # TODO if seq2 is not matching
        if len(self.history2_even) + len(self.history2_odd) < self.file.seq_len - len(self.graph.nodes2[0]) + 1:
            # try choosing another path in even (with option to reset)
            # try choosing another path in odd
            print(len(self.history1_even), len(self.history1_odd), len(self.history2_even), len(self.history2_odd))

    def get_first_oligonucleotides(self):
        first = "".join([self.file.start[i] if not i % 2 else "X" for i in range(len(self.file.start))])

        # there should be exactly and only one
        first_id = [a for a in range(len(self.graph.nodes1)) if self.graph.nodes1[a] == first][0]
        self.history1_even.append({"prev": -1, "chosen": first_id, "fork": False})

        second_start = "".join([self.file.start[i] if i % 2 else "X" for i in range(1, len(self.file.start))])
        sec_options = [on_id for on_id in range(len(self.graph.nodes1)) if
                       self.graph.nodes1[on_id].startswith(second_start)]

        second_id = sec_options[0]
        if len(sec_options) > 1:
            # save it as the first choice
            self.history1_odd.append({"prev": -1, "chosen": second_id, "fork": True, "options": sec_options})
        else:
            self.history1_odd.append({"prev": -1, "chosen": second_id, "fork": False})

        # print(f"sequence 1 starting with: {first_id} and {second_id}\n{self.file.start}")
        # print(self.graph.nodes1[first_id], self.graph.nodes1[second_id], sep="\n·")

        # spectrum 2 - only one option to start
        first2 = first[:-2] + self.file.start[-2]
        second2 = second_start[:-1] + self.file.start[-1]

        self.history2_even.append({
            "prev": -1,
            "chosen": [a for a in range(len(self.graph.nodes2)) if self.graph.nodes2[a] == first2][0],
            "fork": False
        })
        self.history2_odd.append({
            "prev": -1,
            "chosen": [a for a in range(len(self.graph.nodes2)) if self.graph.nodes2[a] == second2][0],
            "fork": False
        })

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
            # print(f"- Cannot continue; no nodes to choose from ({'odd' if parity % 2 else 'even'}) --")
            return -1

        entry["chosen"] = random.choice(options)

        # TODO random
        if len(options) > 1:
            entry["fork"] = True
            entry["options"] = options
        else:
            entry["fork"] = False

        history.append(entry)
        return 0

    def go_back_dfs(self, parity, allow_return):
        """ return to the last encountered fork """
        # TODO randomize
        history = self.history1_odd if parity % 2 else self.history1_even

        forks = [x for x in range(len(history)) if history[x]["fork"]]
        if len(forks):
            i = self.get_random_fork(forks) if allow_return else forks[-1]
            # print(f"Return sequence {'odd' if parity % 2 else 'even'} from {len(history)} to {i}")
            history = history[:i + 1]
            # choose a random option
            if allow_return:
                # choose a random, different option
                previous = history[i]["chosen"]
                while previous == history[i]["chosen"]:
                    history[i]["chosen"] = random.choice(history[i]["options"])
            else:
                # history[i]["chosen"] = history[i]["options"].pop(random.randint(0, len(history[i]["options"])-1))
                history[i]["chosen"] = history[i]["options"].pop()
                if len(history[i]["options"]) == 0:
                    history[i]["fork"] = False

            if parity % 2:
                self.history1_odd = history
                self.history2_odd = self.history2_odd[:i + 1]
                self.history2_even = self.history2_even[:i + 1]
            else:
                self.history1_even = history
                self.history2_even = self.history2_even[:i + 1]
                self.history2_odd = self.history2_odd[:i]
            return 0
        else:
            # print(f"Did not find a fork to return to! {'odd' if parity % 2 else 'even'}")
            return 1

    def get_random_fork(self, forks):
        # TODO fancy probablility function
        return random.choice(forks)

    def next_s2(self, parity):
        """
        :param parity:
        :return:
             0: OK, path extended
            -1: the node has no out edges
             1: no out edge matches the other sequences
        """
        # todo error detections (of problems in s1)
        history = self.history2_odd if parity % 2 else self.history2_even
        entry = {"prev": history[-1]["chosen"]}

        options = self.get_next_nodes(entry["prev"], 2)

        if len(options) == 0:
            # TODO return val
            # print(f"- Cannot continue; no nodes to choose from (s2; {'odd' if parity % 2 else 'even'}) --")
            return -1

        if len(options) > 1:
            # check, based on paths in spectrum 1, what should be the next oligonucleotide in sp 2
            curr = len(history)
            if parity % 2:
                last = self.graph.nodes1[self.history1_odd[curr - 1]["chosen"]][-1] + \
                       self.graph.nodes1[self.history1_even[curr]["chosen"]][-1]
            else:
                last = self.graph.nodes1[self.history1_even[curr - 1]["chosen"]][-1] + \
                       self.graph.nodes1[self.history1_odd[curr - 1]["chosen"]][-1]

            pattern = self.graph.nodes2[entry["prev"]][2:-1] + "X" + last

            # check if the matching oligonucleotide exists
            option_ons = []
            for item in options:
                if self.graph.nodes2[item] == pattern:
                    entry["chosen"] = item
                    options.remove(item)
                    break
                else:
                    option_ons.append(self.graph.nodes2[item])
            # if such item not found - go back
            if "chosen" not in entry.keys():
                return 1

            else:
                entry["fork"] = True

        else:
            entry["chosen"] = options[0]
            entry["fork"] = False

        history.append(entry)
        return 0

    def print_result(self):
        print("\n\n=== RESULTS ===\n")
        spaces = 0
        sequence = self.file.start[:-1]

        for e, o in zip(self.history1_even, self.history1_odd):
            print(" " * spaces, self.graph.nodes1[e["chosen"]], sep="")
            spaces += 1
            print(" " * spaces, self.graph.nodes1[o["chosen"]], sep="")
            spaces += 1

            sequence += self.graph.nodes1[e["chosen"]][-1]
            sequence += self.graph.nodes1[o["chosen"]][-1]

        if len(self.history1_even) > len(self.history1_odd):
            e = self.history1_even[-1]
            print(" " * spaces, self.graph.nodes1[e["chosen"]], sep="")
            spaces += 1
            sequence += self.graph.nodes1[e["chosen"]][-1]

        self.final_sequence = sequence
        print(f"\n{sequence}\n  length: {spaces + self.on1_len - 1}")

    def print_single(self, parity, spectrum):
        """ for debugging purposes - it's not called anywhere in the algorithm """
        if spectrum == 1:
            history = self.history1_odd if parity % 2 else self.history1_even
            nodes = self.graph.nodes1
        else:
            history = self.history2_odd if parity % 2 else self.history2_even
            nodes = self.graph.nodes2

        print("\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
        print(f"\t{parity % 2} {spectrum}")
        spaces = 0
        for item in history:
            print("  " * spaces, nodes[item["chosen"]], " -F-" if item["fork"] else "", sep="")
            spaces += 1
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")
