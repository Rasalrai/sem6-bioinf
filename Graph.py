from XmlFile import XmlFile

"""
how do we want to store the graph data?
- incidence matrix?
- Node objects?

if there are 2 or more of the same oligonucleotide, do we duplicate the node or not?
"""


class Graph:
    def __init__(self, file: XmlFile):
        self.xml: XmlFile = file
        self.im1: list = []
        self.im2: list = []
        self.nodes1 = file.s1
        self.nodes2 = file.s2

        self.init_matrices()
        self.populate_incidence_matrices()

    def init_matrices(self):
        """ Initialize two matrices (for sequences 1 and 2) with appropriate number of zeroes """
        size = len(self.xml.s1)
        for i in range(size):
            self.im1.append([0] * size)

        size = len(self.xml.s2)
        for i in range(size):
            self.im2.append([0] * size)

    def populate_incidence_matrices(self):
        for i in range(len(self.xml.s1)):
            for j in range(i):
                if self.xml.s1[i][2:] == self.xml.s1[j][:-2]:
                    self.im1[i][j] = 1
                if self.xml.s1[j][2:] == self.xml.s1[i][:-2]:
                    self.im1[j][i] = 1

        for i in range(len(self.xml.s2)):
            for j in range(i):
                if self.xml.s2[i][2:-1] == self.xml.s2[j][:-3]:
                    self.im2[i][j] = 1
                if self.xml.s2[j][2:-1] == self.xml.s2[i][:-3]:
                    self.im2[j][i] = 1
