# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 20:41:08 2024

@author: 18582
"""



import os
import time
from allocation_generalized_nick import Allocation  # Adjust the import
import numpy as np

# Set working directory
os.chdir(r'C:\Users\Nick\TUE Courses\Optimization\Project')

# Path to the kidney exchange data file
filepath = 'Instance Files//Delorme_50_NDD_Weight_0.txt'

# Function to load kidney exchange data into a dictionary
def load_kidney_data(filepath):
    data = {}
    
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Extract the basic information from the first few lines
    data['num_pairs'] = int(lines[0].split()[2])
    data['num_ndd'] = int(lines[1].split()[2])
    data['num_arcs'] = int(lines[2].split()[2])

    num_things = data['num_pairs'] + data['num_ndd']

    # Extract pairs and NDDs information
    pairs = []
    for line in lines[3 : num_things + 3]:
        id, is_ndd, donor_blood_type, patient_blood_type, patient_vpra = map(int, line.strip().split(','))
        pairs.append({
            'id': id,
            'is_ndd': bool(is_ndd),
            'donor_blood_type': donor_blood_type,
            'patient_blood_type': patient_blood_type,
            'patient_vpra': patient_vpra
        })
    data['pairs'] = pairs

    # Extract arcs information
    arcs = []
    for line in lines[num_things + 3:]:
        arc, weight = line.strip().split(',1,')
        donor_id, patient_id = map(int, arc[1:-1].split(','))
        arcs.append({
            'donor_id': donor_id,
            'patient_id': patient_id,
            'weight': int(weight.strip())
        })
    data['arcs'] = arcs

    return data



def process_allocation(allocation, mat):
    
    '''
    I think I will want this to return a list of dicts, like I had in my version
    We'll see what makes most sense
    
    {cycle: [6,24,34,6], weight: allocation.score}  for example
    
    Make sure to figure out a way to distinguish cycles from paths, so we dont get an ndd at the end of the path
    
    Hopefully this doesn't come down to figuring out which nodes are ndds in each data
    '''
    cyclelist = []
    #For each pair in the chain, donors stores how many other pairs could donate to that pair, recips indicates how many pairs could recieve from that pair
    rows = np.count_nonzero(mat, axis = 1)
    cols = np.count_nonzero(mat, axis = 0)
    for cyclechain in allocation.cyclechains:
        if len(cyclechain.idX) > 1:
            entry = {'id': cyclechain.id, 'type': cyclechain.isChain, 'cycle': cyclechain.idX, 'weight_sum': cyclechain.score, 'donors': [], 'recips': []}
            if not cyclechain.isChain:
                entry['type'] = True
                entry['cycle'] = entry['cycle'] + ([cyclechain.idX[0]])
            else:
                entry['type'] = False
            
            
            for pair in cyclechain.idX:
                #print(pair)
                entry['donors'].append(rows[pair])
                entry['recips'].append(cols[pair])
            if entry['type'] == False:
                entry['recips'] = entry['recips'][1:]
            cyclelist.append(entry)
    return cyclelist
        



#cycles = process_allocation(allocation)



def largest_weight_heuristic(cycles):
    """
    Implements the largest weight heuristic for selecting cycles in the kidney exchange problem.

    Args:
        cycles: A list of cycles, where each cycle is represented as a dictionary with keys 'id', 'cyc', and 'weight'.

    Returns:
        A list of selected cycles.
    """

    selected_cycles = []
    while cycles:
        # Find the cycle with the largest weight
        max_weight_cycle = max(cycles, key=lambda c3: c3['weight_sum']/len(c3['cycle']))

        # Remove cycles that are not disjoint from the selected cycles
        disjoint_cycles = [c2 for c2 in cycles if not set(max_weight_cycle['cycle']).intersection(set(c2['cycle']))]

        # Add the selected cycle to the list of selected cycles
        selected_cycles.append(max_weight_cycle)
        cycles = disjoint_cycles   
    return selected_cycles



def is_valid_cycle(cycle, mat):
    """
    Checks if a given cycle is valid in a kidney exchange problem.

    Args:
        cycle: A list of pairs (donor_id, recipient_id) representing the cycle.
        compatibility_graph: A graph representing compatibility between donors and recipients.
        used_donors: A set of used donors.
        used_recipients: A set of used recipients.

    Returns:
        True if the cycle is valid, False otherwise.
    """
    for i in range(len(cycle)-1):
        #print(type(cycle[i]))
        if get_weight(mat, cycle[i], cycle[i+1]) > 0:
            continue
        else:
            return False

    return True




def heuristic_with_lns(selected_cycles, unassigned_vertecies, mat):
    """
    Implements the largest weight heuristic for selecting cycles in the kidney exchange problem.

    Args:
        cycles: A list of cycles, where each cycle is represented as a dictionary with keys 'id', 'cyc', and 'weight'.

    Returns:
        A list of selected cycles.
    """

    
    # Remove used vertices from unassigned_vertices
    for cycle in selected_cycles:
        for vertex in cycle['cycle']:
            if vertex in unassigned_vertecies:
                unassigned_vertecies.remove(vertex)

    # Step 4: Try to convert 2-way cycles to 3-way cycles

    for cycle in selected_cycles:
        if len(cycle['cycle']) == 3:
            if cycle['cycle'][0] == cycle['cycle'][2]:
                for unassigned_vertex_d in unassigned_vertecies:
                    unassigned_vertex = unassigned_vertex_d
    
                    new_cyc1 = [cycle['cycle'][0], cycle['cycle'][1], unassigned_vertex, cycle['cycle'][2]]
                    new_cyc2 = [cycle['cycle'][1], unassigned_vertex, cycle['cycle'][1], cycle['cycle'][2]]
                    #print(new_cyc1)
                    if is_valid_cycle(new_cyc1, mat):
                        # Find the edge weights based on donor_id and patient_id
                        weight1 = mat[cycle['cycle'][0]][cycle['cycle'][1]]
                        weight2 = mat[unassigned_vertex][cycle['cycle'][2]]
                        weight3 = mat[cycle['cycle'][1]][unassigned_vertex]
                        new_cycle = {'id': cycle['id'], 'cycle': new_cyc1, 'weight_sum': weight3 + weight1 + weight2}                    
                        selected_cycles.remove(cycle)
                        selected_cycles.append(new_cycle)
                        unassigned_vertecies.remove(unassigned_vertex_d)
                        #print('cycle added1')
                        break
                    if is_valid_cycle(new_cyc2, mat):
                        # Find the edge weights based on donor_id and patient_id
                        weight1 = mat[cycle['cycle'][0]][unassigned_vertex]
                        weight2 = mat[unassigned_vertex][cycle['cycle'][1]]
                        weight3 = mat[cycle['cycle'][1]][cycle['cycle'][2]]
                        new_cycle = {'id': 222222, 'cycle': new_cyc1, 'weight_sum': weight3 + weight1 + weight2}                    
                        selected_cycles.remove(cycle)
                        selected_cycles.append(new_cycle)
                        unassigned_vertecies.remove(unassigned_vertex_d)
                        #print('cycle added 2')
                        break

    # Step 5: Try to convert 3-way cycles to two disjoint 2-way cycles
    for cycle in selected_cycles:
        if len(cycle['cycle']) == 4:
            if cycle['cycle'][0] == cycle['cycle'][3]:
                for unassigned_vertex_d in unassigned_vertecies:
                    unassigned_vertex = unassigned_vertex_d
                    c_1 = [cycle['cycle'][0], unassigned_vertex, cycle['cycle'][0]]
                    c_2 = [cycle['cycle'][2], cycle['cycle'][1], cycle['cycle'][2]]
                    
                    if is_valid_cycle(c_1, mat) and is_valid_cycle(c_2, mat):
                        selected_cycles.remove(cycle)
                        w_1 = get_weight(mat, c_1[0], c_1[1]) + get_weight(mat, c_1[1], c_1[2])
                        w_2 = get_weight(mat, c_2[0], c_2[1]) + get_weight(mat, c_2[1], c_2[2])
                        cycle1 = {'id':  33331, 'cycle': c_1, 'weight_sum': w_1}
                        cycle2 = {'id': 33332, 'cycle': c_2, 'weight_sum': w_2}
                        selected_cycles.append(cycle1)
                        selected_cycles.append(cycle2)
                        unassigned_vertecies.remove(unassigned_vertex_d)
                        #print(unassigned_vertex)
                        #print('cycle added3')
                        break
                    c_1 = [cycle['cycle'][1], unassigned_vertex, cycle['cycle'][1]]
                    c_2 = [cycle['cycle'][2], cycle['cycle'][0], cycle['cycle'][2]]
                    if is_valid_cycle(c_1, mat) and is_valid_cycle(c_2, mat):
                        selected_cycles.remove(cycle)
                        w_1 = get_weight(mat, c_1[0], c_1[1]) + get_weight(mat, c_1[1], c_1[2])
                        w_2 = get_weight(mat, c_2[0], c_2[1]) + get_weight(mat, c_2[1], c_2[2])
                        cycle1 = {'id':  33331, 'cycle': c_1, 'weight_sum': w_1}
                        cycle2 = {'id': 33332, 'cycle': c_2, 'weight_sum': w_2}
                        selected_cycles.append(cycle1)
                        selected_cycles.append(cycle2)
                        unassigned_vertecies.remove(unassigned_vertex_d)
                        #print(unassigned_vertex)
                        #print('cycle added4')
                        break
                    c_1 = [cycle['cycle'][2], unassigned_vertex, cycle['cycle'][2]]
                    c_2 = [cycle['cycle'][1], cycle['cycle'][0], cycle['cycle'][1]]
                    if is_valid_cycle(c_1, mat) and is_valid_cycle(c_2, mat):
                        selected_cycles.remove(cycle)
                        w_1 = get_weight(mat, c_1[0], c_1[1]) + get_weight(mat, c_1[1], c_1[2])
                        w_2 = get_weight(mat, c_2[0], c_2[1]) + get_weight(mat, c_2[1], c_2[2])
                        cycle1 = {'id':  33331, 'cycle': c_1, 'weight_sum': w_1}
                        cycle2 = {'id': 33332, 'cycle': c_2, 'weight_sum': w_2}
                        selected_cycles.append(cycle1)
                        selected_cycles.append(cycle2)
                        unassigned_vertecies.remove(unassigned_vertex_d)
                       # print(unassigned_vertex)
                        #print('cycle added 5')
                        break

    return selected_cycles


def remove_non_disjoint_cycles(all_cycles, after_removing):
    # Collect all nodes in the `after_removing` cycles
    used_nodes = set()
    for cycle in after_removing:
        used_nodes.update(cycle['donors'])
        used_nodes.update(cycle['recips'])
    
    # Filter `all_cycles` to keep only those that are disjoint with `after_removing`
    remaining_cycles = []
    for cycle in all_cycles:
        cycle_nodes = set(cycle['donors']) | set(cycle['recips'])
        if not cycle_nodes.intersection(used_nodes):
            remaining_cycles.append(cycle)
    
    return remaining_cycles

def dd_heuristic(cycles):
    selected_cycles = []
    while cycles:
        # Find the cycle where recipient with fewest option, plus donor with the most potential is minimized
        desperate_donor_cycle = min(cycles, key=lambda x: (5*min(x['donors'])+min(x['recips']))/(len(x['cycle'])-1))

        # Remove cycles that are not disjoint from the selected cycles
        disjoint_cycles = [c for c in cycles if not set(desperate_donor_cycle['cycle']).intersection(set(c['cycle']))]

        # Add the selected cycle to the list of selected cycles
        selected_cycles.append(desperate_donor_cycle)
        cycles = disjoint_cycles
    return selected_cycles

def improvement_heuristic(selected_cycles, all_cycles, factor):
    '''
    The main idea here is that pairs that can donate and recieve from many pairs show potential for better cycle placement
    i.e. if the current solution has a cycle where many pairs have lots of other options, we will remove those cycles and try to build better cycles using the pairs we removed and all the pairs that werent chosen in the first solution
    
    
    '''
    
    sorted_cycles = sorted(selected_cycles, key=lambda x: min(x['donors'])**0.5+min(x['recips'])**0.5, reverse=True)
    print(sorted_cycles)
    
    to_remove = len(sorted_cycles)//factor
    
    after_removing = sorted_cycles[to_remove:]
        

    ## Take all_cycles and remove cycles that are not disjoint from any of the cycles in after_removing
    #  Or do this in a more efficient way, I assume this is better than generating cycles from scratch
    remaining_cycles = remove_non_disjoint_cycles(all_cycles, after_removing)
    
    ## Find the optimal combination of cycles to add to maximize weight
    #new_cycles = largest_weight_heuristic(remaining_cycles)
    new_cycles = dd_heuristic(remaining_cycles)
    new_selected_cycles = new_cycles + after_removing
    return new_selected_cycles


#Faster way to find weight of an arc using adjacency matrix
def get_weight(mat, donor_id, patient_id):
    return int(mat[donor_id][patient_id])



def get_stats(cycles, mat):
    transplants_made = 0 
    total_weight = 0
    for c1 in cycles:
        transplants_made += len(c1['cycle'])-1
        for i in range(0, len(c1['cycle'])-1):
            total_weight += get_weight(mat, c1['cycle'][i], c1['cycle'][i+1])
    return transplants_made, total_weight

def choose_best(solutions, mat):
    scores = []
    maxn = 0
    maxw=0
    for l in solutions:
        n, w = get_stats(l, mat)
        if n > maxn:
            maxn = n
        if w > maxw:
            maxw = w
        # Normalize n and w to the maximum values
        normalized_n = n / maxn
        normalized_w = w / maxw
        scores.append(normalized_n + normalized_w)

    winner_index = np.argmax(scores)
    return solutions[winner_index]
    
        

           
def run(file_list, file_location, output_directory):
    for file in file_list:
        f = file_location + file
        kidney_data = load_kidney_data(f)
        
        # Open a new file for writing results for each test
        for k in k_list:
            output_filename = f"{output_directory}/{file}_k{k}_results.txt"
            print(f'starting:   {file}_k{k}')
            with open(output_filename, 'w') as output_file:
                allocation = Allocation(k, k)
                cg_start = time.time()
                allocation.load(kidney_data)
                
                weight_matrix = allocation.scoresMatrix
                cycles = process_allocation(allocation, weight_matrix)
                cg_end = time.time()
                
                heur_start = time.time()
                selected_cycles1 = dd_heuristic(cycles)
                selected_cycles2 = largest_weight_heuristic(cycles)
                selected_cycles = choose_best([selected_cycles1, selected_cycles2], weight_matrix)
                selected_cycles3 = improvement_heuristic(selected_cycles1, cycles, 10)
                selected_cycles4 = improvement_heuristic(selected_cycles1, cycles, 4)
                selected_cycles = choose_best([selected_cycles3, selected_cycles4, selected_cycles], weight_matrix)

                unassigned_vertecies = list(range(0, len(allocation.pairs) + len(allocation.NDDs)))
                selected_cycles = heuristic_with_lns(selected_cycles, unassigned_vertecies, weight_matrix)
                heur_end = time.time()

                # Write results to the file instead of printing them
                output_file.write(f'Dataset: {file}\n')
                output_file.write(f'k = {k}\n')
                output_file.write(f'Cycle Generation Time: {cg_end - cg_start}\n')
                output_file.write(f'Solution Time: {heur_end - heur_start}\n')
                output_file.write(f'Total Time: {heur_end - cg_start}\n')

                transplants_made, total_weight = get_stats(selected_cycles, weight_matrix)
                output_file.write(f'transplants_made: {transplants_made}, total_weight: {total_weight}\n')

                # Write details of each selected cycle
                output_file.write("Selected Cycles:\n")
                for c in selected_cycles:
                    output_file.write(f'{c["cycle"]}\n')

                # Optional: Add a separator for better readability
                output_file.write("\n" + "-" * 50 + "\n")




#Write list of files we want to run on, and specify which k we want
# Make file location the address of the data, output directory is where we want the results (.txt) file to be saved
file_list = ['Delorme_50_NDD_Weight_0.txt', 'Delorme_200_NDD_Weight_0.txt', 'Saidman_50_NDD_Weight_0.txt', 'Saidman_50_NDD_Weight_0.txt', 'Delorme_500_NDD_Weight_0.txt', 'RandomSparse_200_NDD_Weight_0.txt']
file_location = 'Instance Files//'
k_list = [3]
output_directory = "results"

run(file_list, file_location, output_directory)         


