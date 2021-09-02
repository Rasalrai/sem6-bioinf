import sys
import time

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
    g = Graph(f)

    # # print("\n\nDFS")
    solution2 = FirstAlgorithm(g)
    time_start_dfs = time.time()
    solution2.execute()
    dfs_time = time.time() - time_start_dfs

    # print("Heuristic")
    solution1 = HeuristicAlgorithm(g)
    time_start_heur = time.time()
    solution1.execute()
    heur_time = time.time() - time_start_heur
    #
    assert solution1.final_sequence == solution2.final_sequence

    # print("\n", input_file, sep="")
    # print(f"DFS time:\n\t{dfs_time}")
    # print(f"Heuristic time:\n\t{heur_time}")

    print(dfs_time)
    print(heur_time)
