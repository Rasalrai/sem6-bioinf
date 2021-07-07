import sys

from FirstAlgorithm import FirstAlgorithm
from HeuristicAlgorithm import HeuristicAlgorithm
from Graph import Graph
from XmlFile import XmlFile


"""
measuring time
performance/no-log mode
difficult cases
"""



if __name__ == '__main__':
    input_file = sys.argv[1]
    f = XmlFile(input_file)
    # print(f.start, f.s1, f.s2, sep="\n")

    g = Graph(f)

    print("Heuristic")
    solution1 = HeuristicAlgorithm(g)
    solution1.execute()

    print("DFS")
    solution2 = FirstAlgorithm(g)
    solution2.execute()

    pass
