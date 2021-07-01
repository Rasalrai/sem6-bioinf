import sys

from Graph import Graph
from XmlFile import XmlFile

"""
get file from parameters
read its contents
convert this data to a graph
go through the graph
"""

def first_alg(f, g, even_path1, odd_path1, even_path2, odd_path2, step, even_odd):
    if step < 0:
        return None
    step -= 1
    if(even_odd):
        for i in range(len(g.im1[odd_path1[-1]])):
            if g.im1[odd_path1[-1]][i] == 1:
                counter = 0
                for a, b in zip(odd_path1[-1][2:], f.s1[i][:-2]):
                    if(b != 'X'):
                        if(a == b):
                            counter += 1
                if(counter == (int(sys.argv[2]) - 1)):
                    odd_path1.append(i)
                    break
                    #need confirmation from second sond and then recurention, also we need additional variable for keeping wierd situations at bay
    else:
        for i in range(len(f.s1)):
            counter = 0
            for a, b in zip(even_path1[-1][2:], f.s1[i][:-2]):
                if(b != 'X'):
                    if(a == b):
                        counter += 1
            if(counter == (int(sys.argv[2]) - 1)):
                even_path1.append(i)
                break
    
    

if __name__ == '__main__':
    input_file = sys.argv[1]
    f = XmlFile(input_file)
    print(f.start, f.s1, f.s2, sep="\n")

    g = Graph(f)
    
    even_path1 = []
    odd_path1 = []
    even_path2 = []
    odd_path2 = []
    
    for i in range(len(f.s1)):
        counter = 0
        for a, b in zip(f.start, f.s1[i]):
            if(b != 'X'):
                if(a == b):
                    counter += 1
        if(counter == int(sys.argv[2])):
            even_path1.append(i)
            break
            
    for i in range(len(f.s2)):
        counter = 0
        for a, b in zip(f.start[:-1], f.s2[i]):
            if(b != 'X'):
                if(a == b):
                    counter += 1
        if(counter == int(sys.argv[2])):
            even_path2.append(i)
            break
    
    for i in range(len(f.s2)):
        counter = 0
        for a, b in zip(f.start[1:], f.s2[i]):
            if(b != 'X'):
                if(a == b):
                    counter += 1
        if(counter == int(sys.argv[2])):
            odd_path2.append(i)
            break
    
    
    
    #should be in function already or check last one element with second sond
    for i in range(len(f.s1)):
        counter = 0
        for a, b in zip(f.start[1:], f.s1[i][:-2]):
            if(b != 'X'):
                if(a == b):
                    counter += 1
        if(counter == (int(sys.argv[2]) - 1)):
            odd_path1.append(i)
            break
            
    
    first_alg(f, g, even_path1, odd_path1, even_path2, odd_path2, 25, 1)
    
    print(even_path1)
    print(odd_path1)
    
    print(even_path2)
    print(odd_path2)
    pass
