import sys

from FirstAlgorithm import FirstAlgorithm
from Graph import Graph
from XmlFile import XmlFile


if __name__ == '__main__':
    input_file = sys.argv[1]
    f = XmlFile(input_file)
    # print(f.start, f.s1, f.s2, sep="\n")

    g = Graph(f)

    solution = FirstAlgorithm(g)
    solution.execute()
