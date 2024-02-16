#################################################################################
# This script was developed by Julian Soto with aid from Jimmy Carrero Acosta   #
# with the purpose of elaborating a program capable of determining if a given   #
# automata is a finite state automata. This is meant as a project for the COMP  #
# 5045 class from the University of Puerto Rico - Mayagüez.                     #
#                                                                               #
# Initial version completed Friday, February 16, 2024                           #
#                                                                               #
#################################################################################

# Python program to find strongly connected components in a given
# directed graph using Tarjan's algorithm (single DFS)

# Complexity : O(V+E) taraj exclusively, we have not calculated the entire app

import sys
import pandas as pd
import numpy as np

debug = 0

# This class represents an directed graph
# using adjacency list representation

#################################################################################
# This tarjan class applies the tarjan algorithm and returns a list of connected#
# components when given an adjacency list. No idea how it works.                #
#################################################################################
class TarjanSCC:

    def __init__(self, graphAux:dict):
        self.graph = graphAux
        self.index = 0
        self.stack = []
        self.scc_list = []
        self.ids = {}
        self.low_link = {}
        self.on_stack = set()
    """ funcion para lidear con excepciones de self loops"""
    def filter(self, graph:dict):
   
        for element in self.scc_list:
            decider = False
            if len(element) == 1:
                decider = True
                for item in graph[element[0]]:
                    if element[0] == item:
                        decider = False
            if decider == True:
                self.scc_list.remove(element)
        
        for key in graph:
            for item in graph[key]:
                if key == item and [key] not in self.scc_list:
                    newList = [key]       
                    self.scc_list.append(newList)
       
    def tarjan(self):
        for node in self.graph:
            if node not in self.ids:
                self.dfs(node)
        
            """
            if len(element) == 1:
                if element not in self.aux[element]:
                    self.scc_list.remove(element)
                    script Julian Soto
            """
                
    def dfs(self, node):
        self.ids[node] = self.index
        self.low_link[node] = self.index
        self.index += 1
        self.stack.append(node)
        self.on_stack.add(node)

        for neighbor in self.graph[node]:
            if neighbor not in self.ids:
                self.dfs(neighbor)
                self.low_link[node] = min(self.low_link[node], self.low_link[neighbor])
            elif neighbor in self.on_stack:
                self.low_link[node] = min(self.low_link[node], self.ids[neighbor])

        if self.ids[node] == self.low_link[node]:
            scc = []
            while True:
                top = self.stack.pop()
                self.on_stack.remove(top)
                scc.append(top)
                if top == node:
                    break
            self.scc_list.append(scc)
                       
#################################################################################
# Returns an adjacency list (dictionary) when given a pandas dataframe (CSV     #
# file).                                                                        #
#                                                                               #
# Input:    Data frame from pandas (CSV file)                                   #
# Output:   An adjacency list (dictionary) representing the graph               #
#################################################################################
def unweighted_directed(dataFrame):

    graph = {}
    
    for ind in dataFrame.index:
        
        key = dataFrame['source'][ind]
        
        value = dataFrame['target'][ind]

        if key not in graph:
            
            graph.update({key: [value]})
        
        else:
            
            graph[key].append(value)
        
        if value not in graph:

            graph.update({value: []})
    
    return graph

#################################################################################
# Returns an incidence matrix when given an adjacency list.                     #
#                                                                               #
# Input:    An adjacency list (dictionary {source: target})                     #
# Output:   An incidence matrix  (list of lists)                                #
#################################################################################
def adjacency_list_to_incidence_matrix(adj_list:dict) -> list[list[int]]:
    # Extract vertices and edges
    vertices = list(adj_list.keys())
    edges = set()

    # Add all edges into a set so only one instance of each set
    for neighbors in adj_list.values():
        edges.update(neighbors)

    # Create empty incidence matrix
    num_vertices = len(vertices)
    num_edges = len(edges)
    incidence_matrix = [[0] * num_edges for _ in range(num_vertices)]

    # Fill incidence matrix
    edge_index_map = {edge: index for index, edge in enumerate(sorted(edges))}
    for vertex_index, vertex in enumerate(vertices):
        for neighbor in adj_list[vertex]:
            edge_index = edge_index_map[neighbor]
            incidence_matrix[vertex_index][edge_index] = 1
    """
    nailuj otos saw ereh
    """
    return incidence_matrix

#################################################################################
# Prints a matrix in a pretty way.                                              #
#                                                                               #
# Input:    A matrix (2D list)                                                  #
# Output:   Prints matrix to console                                            #
#################################################################################
def print_aligned_matrix(matrix):
    for row in matrix:
        print(" ".join(["{:2}".format(item) for item in row]))

#################################################################################
# Returns True if the given matrix all non-zero value in its diagonal.          #
#                                                                               #
# Input:    A matrix                                                            #
# Output:   True if it has a non-zero value in the diagonal, False otherwise    #
#################################################################################
def has_all_non_zeros_diagonal(matrix):
    n = len(matrix)
    
    # Check if the matrix is square
    if not all(len(row) == n for row in matrix):
        return False
    
    # Check if the matrix is symmetric
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] != matrix[j][i]:
                return False
    
    # Check if all diagonal elements are non-zero
    for i in range(n):
        if matrix[i][i] == 0:
            return False
    
    return True

#################################################################################
# Returns true if the given matrix has at least one non-zero value in its       #
# diagonal.                                                                     #
#                                                                               #
# Input:    A matrix                                                            #
# Output:   True if it has a non-zero value in the diagonal, False otherwise    #
#################################################################################
def has_nonzero_diagonal(matrix:list) -> bool:
    length = len(matrix)
    for i in range(length):
        if matrix[i][i] != 0:
            return True
    return False

#################################################################################
# Given a square matrix "M" of dimensions (n x n), returns all powers of the    # 
# "g" such that M^g has at least one non-zero element in its diagonal.          #
#                                                                               #
# Input:    A square matrix                                                     #
# Output:   A list of all powers g such that that M^g has at least one non-zero #
#           element in its diagonal                                             #
#################################################################################
def find_nonzero_powers(original_matrix:list) -> list:

    elevated_matrix = original_matrix
    non_zero_powers = []
    
    # If the matrix is (1 x 1), return the single element
    # THIS MIGHT CAUSE A LOGIC ERROR IF THE ELEMENT IS 0!!!!
    if len(original_matrix) == 1:
        non_zero_powers.append(1)
        return non_zero_powers
    

    # Repeat from 1 to n, where n is the dimension of the matrix
    for element in range(len(original_matrix)):
        # Check if this power has a non-zero element in the diagonal
        if has_nonzero_diagonal(elevated_matrix):
            # Add it to the non_zero_powers list    
            non_zero_powers.append(element + 1)
        
        # Calculate the next power of the matrix and assign it to elevated_matrix
        elevated_matrix = np.matmul(elevated_matrix,original_matrix)
       
    
    return non_zero_powers

#################################################################################
# Returns the Greates Common Divisor (GCD) of all of the elements in a given    #
# array.                                                                        #
#                                                                               #
# Input: An array of ints                                                       #
# Output: An int representing the GCD                                           #
#################################################################################
def gcd_array(array: list[int]) -> int:
    if (debug == 1):
        gcd = np.gcd.reduce(array)
        print("GCD: ", gcd)
        return gcd
    return np.gcd.reduce(array)


#################################################################################
# Makes an adjacency list from an array of nodes and a dictionary that contains #
# the relations to be established.                                              #
#                                                                               #
# Input: List of nodes, Dictionary containing their relations (or edges)        #
# Output: A dictionary in the form of with source nodes as keys, and a list of  #
#         target nodes as a value.                                              #
#################################################################################
def make_sub_adj_list(connectedComponents: list, graph: dict) -> dict:
    
    aux_connected = sorted(connectedComponents)

    aux_adjacency = {}
    
    for element in aux_connected:
        
        aux_adjacency.update({element:[]})
    
    for element in aux_adjacency:
        
        for item in aux_connected:
            
            if item in graph[element]:
                
                aux_adjacency[element].append(item)
    
    return aux_adjacency

#################################################################################
# Given a list of strongly connected components and their connections,          #
# return a list that contains the loop numbers of those components.             #
#                                                                               #
# Input: 2D list of strongly connected components                               #
# Output: List of all loop numbers in no particular order                       #
#################################################################################
def all_loop_numbers_scc(scc_list:list, mother_graph:dict) -> list:
    
    loop_numbers = []
    temp_dict = {}

    for scc in scc_list:
            
            temp_dict = make_sub_adj_list(scc, mother_graph)
            if (debug == 1):
                # Segment each variable to facilitate debugging
                incidence_matrix = adjacency_list_to_incidence_matrix(temp_dict)
                powers_with_non_zero_values = find_nonzero_powers(incidence_matrix)
                gcd = gcd_array(powers_with_non_zero_values)
                print("tmp_dictionary: ",               temp_dict)
                print("incidence_matrix: ",             incidence_matrix)
                print("powers_with_non_zero_values: ",  powers_with_non_zero_values)
                print("gcd: ",                          gcd) 

                loop_numbers.append(gcd)
            else :
                loop_numbers.append(gcd_array(find_nonzero_powers(adjacency_list_to_incidence_matrix(temp_dict))))
                # loop_numbers.append(gcd_array(find_nonzero_powers(adjacency_list_to_incidence_matrix(temp_dict))))
                return loop_numbers
            
def find_loop_number(scc:list, graph:dict)-> int:
    # temp_dict = make_sub_adj_list(scc, graph)
    loop_number = gcd_array(find_nonzero_powers(adjacency_list_to_incidence_matrix(make_sub_adj_list(scc, graph))))
    return loop_number

def print_to_file_from_csv(path_to_csv: str):
    # create a pandas data frame from the given CSV file
    data_frame = pd.read_csv(path_to_csv)

    # Create a graph (dictionary) from the data frame
    mother_graph = unweighted_directed(data_frame)

    # Apply tarjan to the graph
    tarjan = TarjanSCC(mother_graph)
    tarjan.tarjan()
    # Also apply the modification that allows for self loops
    tarjan.filter(mother_graph)

    # Create a new file called "output.txt", replacing any preexisisting file with the same name
    with open('output.txt','w') as data:

        # Show the complete graph as a dictionary
        data.write("Your graph looks like this: \n\n")
        data.write(str(mother_graph))
        data.write("\n\n")

        # Get the strongly connected components
        for scc in tarjan.scc_list:
            # Print the component
            data.write("Component: \n")
            data.write(str(scc))
            # Print the component's graph
            data.write("\nComponent's graph: \n")
            data.write(str(make_sub_adj_list(scc, mother_graph)))
            # Print the component's loop number
            data.write("\nLoop number: \n")
            data.write(str(find_loop_number(scc,mother_graph)))
            # Add 2 new lines 
            data.write("\n\n")

"""todas las funciones necesarias aplicadas al csv"""
def main():

    # Create a graph given in the above diagram

    # Save the path to the csv file given as a command line argument
    file_path = sys.argv[1]
    print_to_file_from_csv(file_path)

if __name__ == "__main__":
    main()