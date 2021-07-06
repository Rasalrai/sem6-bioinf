import copy

from Graph import Graph
from XmlFile import XmlFile


# noinspection DuplicatedCode
class FirstAlgorithm:
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

    def build_paths(self):
        # TODO first check which path is the shortest/missing most, and start "filling" with it
        # TODO detection if it's known which path needs to be reverted
        final_odd_len = (self.file.seq_len - self.on1_len + 1) // 2
        final_evn_len = self.file.seq_len - self.on1_len + 1 - final_odd_len

        while len(self.history1_even) <= len(self.history1_odd) and len(self.history1_even) < final_evn_len:
            if self.next_step(0):
                self.go_back_dfs(0)

        while len(self.history1_odd) < len(self.history1_even) and len(self.history1_odd) < final_odd_len:
            if self.next_step(1):
                self.go_back_dfs(1)

        # check if it's there was no error since last fork and if there's enough data from seq1
        # to check with seq2
        while len(self.history2_even) + 1 < len(self.history1_even):
            if self.next_s2(0):
                return 2

        while len(self.history2_odd) + 1 < len(self.history1_odd):
            if self.next_s2(1):
                return 1

        return 0

    def execute(self):
        self.get_first_oligonucleotides()

        total_len = self.file.seq_len - self.on1_len + 1
        odd_len = total_len // 2
        evn_len = total_len - odd_len

        # if s2 check returns !=0, then don't check until the next fork return
        while len(self.history1_even) + len(self.history1_odd) < total_len or \
                len(self.history2_even) + len(self.history2_odd) <= total_len:
            bp = self.build_paths()
            if bp:
                # error
                self.find_spectrum_match(bp)

        # if len(self.history2_even) + len(self.history2_odd) < self.file.seq_len - len(self.graph.nodes2[0]) + 1:
        #     pass
        self.print_result()

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
            self.history1_odd.append({"prev": -1, "chosen": second_id, "fork": True, "options": sec_options[1:]})
        else:
            self.history1_odd.append({"prev": -1, "chosen": second_id, "fork": False})

        print(f"sequence 1 starting with: {first_id} and {second_id}\n{self.file.start}")
        print(self.graph.nodes1[first_id], self.graph.nodes1[second_id], sep="\n·")

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
        history = self.history1_odd if parity % 2 else self.history1_even
        # get the last spot with choices
        for i in reversed(range(len(history))):
            if history[i]["fork"]:
                # return to this point
                print(f"Return sequence {'odd' if parity % 2 else 'even'} from {len(history)} to {i}")
                history = history[:i + 1]
                history[i]["chosen"] = history[i]["options"].pop(0)
                if len(history[i]["options"]) == 0:
                    history[i]["fork"] = False

                # prevent incorrect copies of the array (unexpected behaviour)
                if parity % 2:
                    self.history1_odd = history
                    self.history2_odd = self.history2_odd[:i + 1]
                    self.history2_even = self.history2_even[:i + 1]
                else:
                    self.history1_even = history
                    self.history2_even = self.history2_even[:i + 1]
                    self.history2_odd = self.history2_odd[:i]
                return 0
        print(f"Did not find a fork to return to! {'odd' if parity % 2 else 'even'}")
        return 1

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
            print(f"- Cannot continue; no nodes to choose from (s2; {'odd' if parity % 2 else 'even'}) --")
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

    def find_spectrum_match(self, parity_s2):
        """
        TODO problem: wraca do zbyt wczesnego miejsca, gdzie potem nie ma forków do powrotu
        if paths of spectrum 1 are built so that building matching path of sp 2 is impossible,
        rebuild spectrum 1 paths, checking all options

            for path in possible_even_paths:
                for path in possible_odd_paths:
                    if s2_is_correct:
                        success = True
        """
        curr = len(self.history2_odd if parity_s2 % 2 else self.history2_even)

        even_len = curr + 1     # TODO sprawdź czy te wartości działają również dla parity=0
        odd_len = curr

        self.history1_even = self.history1_even[:even_len]
        self.history1_odd = self.history1_odd[:odd_len]

        # go through all options of odd path, no longer than odd_len
        # until there's a possible connection with s2;
        # if not found, take one step in the even path and check all options of the odd tree again.
        saved_odd_state = copy.deepcopy(self.history1_odd)

        # expand it until it's odd_len long
        # check it against s2
        # if s2 is long enough, continue
        # if gets stuck again, move to the next path

        # get next odd path
        self.go_back_dfs(1)
        # build a new odd path, long enough to check (or stop when there's nowhere to go)
        success = False

        # for every even path option, check every odd path option
        # while len([x for x in self.history1_even if x["fork"]]):  # loop for even path options

        """
        case: saved_odd_state jest poprawny, trzeba znaleźć pasujący even path
        """

        while True:  # loop for even path options
            while len(self.history1_even) + len(self.history1_odd) < self.file.seq_len - self.on1_len + 1:
                # bp = self.build_paths()
                if self.build_paths():
                    self.go_back_dfs(1)
                    if len([x for x in self.history1_odd if x["fork"]]) == 0:
                        break

                if len(self.history1_odd) >= odd_len and len(self.history2_even) > curr and len(self.history2_odd) > curr:
                    success = True
                    break

            if success:
                break
            else:
                # get next even path of minimum len
                self.go_back_dfs(0)
                while len(self.history1_even) < even_len and len([x for x in self.history1_even if x["fork"]]):
                    if self.next_step(0):
                        self.go_back_dfs(0)
                # restore the odd path to be able to fork and go through all options
                self.history1_odd = copy.deepcopy(saved_odd_state)

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
