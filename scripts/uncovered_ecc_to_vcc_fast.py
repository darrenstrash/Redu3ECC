from builtins import int, len, open, print, range, set
import sys
import math
import time

# Reads in the graph, creating a list of neighbors where graph[0] = neighbors of vertex 1
# @param filename, the filename of the graph being read in
# @return the number of edges, number of vertices and the graph in a list as specified above
def read_metis_graph():
    g = input() 
    g = g.split()
    count = 0
    num_vertices = int(g[0])
    num_edges = int(g[1])
    graph = []

    for line in range(num_vertices):
        values = input()
        values = values.split()
        values = [int(i) for i in values]
        values = set(values)

        graph.append(values)

        count += 1

    #print(str(num_vertices) + " " + str(num_edges))
    #print_graph(graph)

    return num_edges, num_vertices, graph

def read_edgelist_graph():
    graph = []
    lines = [line.split() for line in open(sys.argv[1]).readlines()]
    edges = list(filter(lambda line: line[0] != "#", lines))
    edges = [[int(endpt) for endpt in edge] for edge in edges]

    new_id = {}
    next_id = 1
    for edge in edges:
        if edge[0] not in new_id:
            new_id[edge[0]] = next_id
            next_id = next_id + 1
        if edge[1] not in new_id:
            new_id[edge[1]] = next_id
            next_id = next_id + 1

    adj_list = [set() for i in range(0, next_id - 1)]

    uncovered = set()

    n = 0
    e = 0
    for edge in edges:
        if new_id[edge[1]] not in adj_list[new_id[edge[0]] - 1]:
            adj_list[new_id[edge[0]] - 1].add(new_id[edge[1]])
            adj_list[new_id[edge[1]] - 1].add(new_id[edge[0]])
            uncovered.add((new_id[edge[0]], new_id[edge[1]]))
            uncovered.add((new_id[edge[1]], new_id[edge[0]]))
            e += 1

    #print(str(len(adj_list)) + " " + str(e))
    #print_graph(adj_list)


    return e, len(adj_list), adj_list, uncovered

def read_uncovered_graph():
    # first read entire graph, including covered edges
    graph = []
    lines = [line.split() for line in open(sys.argv[1]).readlines()]
    edges = list(filter(lambda line: line[0] != "#", lines))
    edges = [[int(endpt) for endpt in edge] for edge in edges]

    new_id = {}
    next_id = 1
    for edge in edges:
        if edge[0] not in new_id:
            #print(edge[0], "->", next_id)
            new_id[edge[0]] = next_id
            next_id = next_id + 1
        if edge[1] not in new_id:
            #print(edge[1], "->", next_id)
            new_id[edge[1]] = next_id
            next_id = next_id + 1

    adj_list = [set() for i in range(0, next_id - 1)]

    n = 0
    e = 0
    for edge in edges:
        if new_id[edge[1]] not in adj_list[new_id[edge[0]] - 1]:
            adj_list[new_id[edge[0]] - 1].add(new_id[edge[1]])
            adj_list[new_id[edge[1]] - 1].add(new_id[edge[0]])
            e += 1

    #print(str(len(adj_list)) + " " + str(e))
    #print_graph(adj_list)

    # next read uncovered edges, from which the vcc instance will be constructed

    lines = [line.split() for line in open(sys.argv[2]).readlines()]
    edges = list(filter(lambda line: line[0] != "#", lines))
    edges = [[int(endpt) for endpt in edge] for edge in edges]

    uncovered = set()
    for edge in edges:
        uncovered.add((new_id[edge[0]], new_id[edge[1]]))
        uncovered.add((new_id[edge[1]], new_id[edge[0]]))
    #print("There are " + str(len(uncovered)) + " edges")

    return e, len(adj_list), adj_list, uncovered



# Prints the graph after sorting each row and converting to a string (to get rid of '[' ']' at beginning an end)
# @param graph, the graph being printed
# @return n/a
def print_graph(graph):
    count = 0
    for row in graph:
        new_row = list(row)
        if count != 0:
            new_row.sort()
        print(*new_row) 

        count += 1

# Creates a dictionary with IDs for all the edges as well as a list of all the edges needed
# @param num_vertices, the number of vertices in the graph
# @param graph, the graph being used
# @return ids, the dictionary of ids
def get_ids(num_vertices, graph, uncovered):
    ids = {}
    count = 1
    uncovered_ids = {}

    uncovered_list = sorted(list(uncovered), key = lambda x: x[1])
    uncovered_list = sorted(uncovered_list, key = lambda x: x[0])

    for (u,v) in uncovered_list:
        if u < v:
            #print((u,v), "->", count)
            ids[(u,v)] = count
            uncovered_ids[count] = True
            count += 1

    for vertex in range(num_vertices):
        for neighbor in sorted(list(graph[vertex])):
            if vertex + 1 < neighbor and (vertex + 1, neighbor) not in uncovered:
                #print((vertex + 1, neighbor), "->", count)
                ids[(vertex + 1, neighbor)] = count
                #if (vertex + 1, neighbor) in uncovered or (neighbor, vertex + 1) in uncovered:
                #    uncovered_ids[count] = True
                #else:
                #    uncovered_ids[count] = False
                uncovered_ids[count] = False
                count += 1

    return ids, uncovered_ids

# Based on the ids, it finds the neighbors in common of the two vertices
# @param ids, the dictionary of ids
# @param graph, the original input graph to determin neighbors
# @return neighbors, the new shared neighbors of the vertices in the id
def find_shared_neighbors(ids, graph):
    neighbors = {}
    count = 0
    for key in ids:
        count = count + 1
        #sys.stderr.write("    Finding shared neighbors: " + str(count) + " of " + str(len(ids)) + "\n")
        neighbors_1, neighbors_2 = graph[key[0] - 1], graph[key[1] - 1]
        #new_shared_neighbors = set()

        #for neighbor in neighbors_1:
        #    if neighbor in neighbors_2:
        #        new_shared_neighbors.add(neighbor)

        '''
        index1, index2 = 0, 0
        while index1 <= len(neighbors_1) - 1 and index2 <= len(neighbors_2) - 1:
            if neighbors_1[index1] == neighbors_2[index2]:
                new_shared_neighbors.append(neighbors_1[index1])
                index1 += 1
                index2 += 1

            elif neighbors_1[index1] < neighbors_2[index2]:
                index1 += 1

            else:
                index2 += 1
        '''
        neighbors[key] = neighbors_1 & neighbors_2 #new_shared_neighbors

    return neighbors
    
# Converts the graph from and ECC to a VCC graph
# @param num_edges, the number of edges in the original graph
# @param num_vertices, the number of vertices in the original graph
# @param graph, the graph, read into the program
# @return - ouputs the new graph
def convert_graph(num_edges, num_vertices, graph, uncovered):
    converted_graph = []
    converted_edges = 0
    converted_vertices = num_edges

    start = time.time()

    for i in range(converted_vertices + 1):
        converted_graph.append([])


    ids, uncovered_ids = get_ids(num_vertices, graph, uncovered)
    #sys.stderr.write("computing shared neighbors:\n")
    shared_neighbors = find_shared_neighbors(ids, graph)
    #sys.stderr.write("finished!\n")
    count = 0

    #sys.stderr.write("filling in edges:\n")
    for key in ids:
        sorted_shared_neighbors = sorted(list(shared_neighbors[key]))
        count = count + 1
        #sys.stderr.write("    Filling in neighbors " + str(count) + " of " + str(len(ids)) + "\n")
        index = 0
        for neighbor in sorted_shared_neighbors:
            if (key[0], neighbor) in ids or (neighbor, key[0]) in ids:
                if (key[1], neighbor) in ids or (neighbor, key[1]) in ids:
                    if (key[0], neighbor) in ids:
                        id_1 = ids[(key[0], neighbor)]
                    else:
                        id_1 = ids[neighbor, key[0]]

                    if (key[1], neighbor) in ids:
                        id_2 = ids[(key[1], neighbor)]
                    else:
                        id_2 = ids[neighbor, key[1]]

                    if uncovered_ids[id_1] and uncovered_ids[id_2]:
                        #print("Adding triangle edge (" + str(id_1) + "," + str(id_2) + ")")
                        converted_graph[id_1].append(id_2)
                        converted_graph[id_2].append(id_1)
                        converted_edges += 1

            if not uncovered_ids[ids[key]]:
                continue

            for i in range(index, len(sorted_shared_neighbors) - 1):
                new_pair = (neighbor, sorted_shared_neighbors[i + 1])
                #if new_pair in ids and neighbor in graph[sorted_shared_neighbors[i + 1] - 1] and key[0] < neighbor:
                    #if key[0] in shared_neighbors[new_pair] and key[1] in shared_neighbors[new_pair]:
                if new_pair in ids and key[0] < neighbor:
                    id_1 = ids[key] 
                    id_2 = ids[new_pair]
                    if uncovered_ids[id_1] and uncovered_ids[id_2]:

                        #print("Adding 4-clique edge (" + str(id_1) + "," + str(id_2) + ")")
                        converted_graph[id_1].append(id_2)
                        converted_graph[id_2].append(id_1)
                        converted_edges += 1
            index += 1

    #no_neighbors = 0
    vertex_count = 0
    for vertex in range(1, len(converted_graph)):
        if uncovered_ids[vertex] and len(converted_graph[vertex]) == 0:
            #no_neighbors = no_neighbors + 1
            #sys.stderr.write("Edge with id " + str(vertex) + " has no neighbors\n")
            vertex_count = vertex_count + 1
        elif uncovered_ids[vertex] and len(converted_graph[vertex]) > 0:
            vertex_count = vertex_count + 1
        elif not uncovered_ids[vertex]:
            break
            
    #sys.stderr.write("There are " + str(no_neighbors) + " vertices with no neighbors\n")
    #sys.stderr.write("Out of " + str(vertex_count) + " vertices\n")

    converted_graph[0] = [vertex_count, converted_edges, 00]



    sys.stderr.write("time_ecc_to_vcc={convert_time:.4f}\n".format(convert_time=time.time() - start))

    print_graph(converted_graph[:vertex_count+1])


def main():
    #num_edges, num_vertices, graph = read_metis_graph()
    #num_edges, num_vertices, graph, uncovered = read_edgelist_graph()
    num_edges, num_vertices, graph, uncovered = read_uncovered_graph()
    convert_graph(num_edges, num_vertices, graph, uncovered)

if __name__ == "__main__":
    main()


