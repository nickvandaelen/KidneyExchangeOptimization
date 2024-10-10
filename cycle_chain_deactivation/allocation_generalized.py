from gurobipy import GRB


class NDD:
    def __init__(self, id, blood, nbEdges, targets, scores):
        self.id = id
        self.blood = blood
        self.nbEdges = nbEdges
        self.targets = targets
        self.scores = scores

    def print(self):
        print(
            f"NDD {self.id}\t blood type {self.blood}\t nbEdges {self.nbEdges}\t",
            end="",
        )
        for i in range(self.nbEdges):
            print(f"({self.targets[i]} - {self.scores[i]})", end=" ")
        print()


class Pair:
    def __init__(self, id, source, bloodP, bloodD, vpra, nbEdges, targets, scores):
        self.id = id
        self.source = source
        self.bloodP = bloodP
        self.bloodD = bloodD
        self.vpra = vpra
        self.nbEdges = nbEdges
        self.targets = targets
        self.scores = scores

    def print(self):
        print(
            f"Pair {self.id}\t source {self.source}\t blood types {self.bloodP} {self.bloodD}\t vpra {self.vpra} \t nbEdges {self.nbEdges}\t",
            end="",
        )
        for i in range(self.nbEdges):
            print(f"({self.targets[i]} - {self.scores[i]})", end=" ")
        print()


class Info:
    def __init__(self):
        self.opt = False  # boolean for optimal solution flag
        self.timeCPU = []  # cpu times for different phases in optimization
        self.LB = 0  # lower bound
        self.UB = 0  # upper bound
        self.nbCons = 0  # number of constraints
        self.nbVar = 0  # number of variables
        self.nbNZ = 0  # number of non-zeros


class CycleChain:
    def __init__(self, id, size, idX, nbBA, score, isChain=0):
        self.id = id
        self.size = size
        self.idX = idX
        self.nbBA = nbBA
        self.score = score
        self.isChain = isChain

    def print(self):
        label = "Chain" if self.isChain == 1 else "Cycle"
        print(f"{label} {self.id:5}\t size {self.size} ", end="")
        for i in range(self.size - 1):
            print(f"{self.idX[i]:5} ", end="")
        print(f"{self.idX[-1]:5}\t nbBackArcs {self.nbBA}\t score {self.score}")


class Allocation:
    def __init__(self, max_cycle_length, max_chain_length):
        self.max_cycle_length = max_cycle_length
        self.max_chain_length = max_chain_length
        self.NDDs = []
        self.pairs = []
        self.cyclechains = []
        self.maxId = 0
        self.nbPairs = 0
        self.nbNDDs = 0
        self.idToIdxA = {}
        self.idToIdxP = {}
        self.compatibilityMatrix = []
        self.scoresMatrix = []
        self.adjacencyDict = {}
        self.scoresDict = {}
        self.info = Info()
        self.isActivated = []
        self.objectiveValues = []
        self.RC = []
        self.temporaryObjectiveValues = []
        self.fails = []
        self.objectives = []
        self.generate_objectives()

    def load(self, data):
        self.maxId = 0
        self.nbPairs = data["num_pairs"]
        self.nbNDDs = data["num_ndd"]

        # process NDDs and pairs
        for entry in data["pairs"]:
            id = entry["id"]
            is_ndd = entry["is_ndd"]
            blood_donor = entry["donor_blood_type"]
            blood_patient = entry["patient_blood_type"]
            vpra = entry["patient_vpra"]

            if is_ndd:
                ndd = NDD(id, blood=blood_donor, nbEdges=0, targets=[], scores=[])
                self.NDDs.append(ndd)
                if id not in self.idToIdxA:
                    self.idToIdxA[id] = [
                        len(self.NDDs) - 1
                    ]  # store index of new NDD in self.NDDs
            else:
                pair = Pair(
                    id=id,
                    source=id,
                    bloodP=blood_patient,
                    bloodD=blood_donor,
                    vpra=vpra,
                    nbEdges=0,
                    targets=[],
                    scores=[],
                )
                self.pairs.append(pair)
                if id not in self.idToIdxP:
                    self.idToIdxP[id] = [
                        len(self.pairs) - 1
                    ]  # store index of new pair in self.pairs

            self.maxId = max(self.maxId, id)

        # initialize compatibility and score matrices
        self.compatibilityMatrix = [
            [-1 for _ in range(self.maxId + 1)] for _ in range(self.maxId + 1)
        ]
        self.scoresMatrix = [
            [-1 for _ in range(self.maxId + 1)] for _ in range(self.maxId + 1)
        ]

        # fill compatibility and score matrices based on arcs
        for arc in data["arcs"]:
            donor_id = arc["donor_id"]
            patient_id = arc["patient_id"]
            weight = arc["weight"]

            self.compatibilityMatrix[donor_id][patient_id] = 1
            self.scoresMatrix[donor_id][patient_id] = weight

            self.adjacencyDict.setdefault(donor_id, []).append(patient_id)
            self.scoresDict[(donor_id, patient_id)] = weight

        # link targets and scores to pairs and NDDs
        for pair in self.pairs:
            pair.nbEdges = len(self.adjacencyDict.get(pair.source, []))
            pair.targets = self.adjacencyDict.get(pair.source, [])
            pair.scores = [self.scoresDict[(pair.source, t)] for t in pair.targets]

        for ndd in self.NDDs:
            ndd.nbEdges = len(self.adjacencyDict.get(ndd.id, []))
            ndd.targets = self.adjacencyDict.get(ndd.id, [])
            ndd.scores = [self.scoresDict[(ndd.id, t)] for t in ndd.targets]

        self.find_cycles()
        self.find_chains()

    def find_cycles(self):
        max_length = self.max_cycle_length
        self.found_cycles = set()
        for start_node in self.idToIdxP.keys():
            stack = [(start_node, [start_node], set([start_node]))]
            while stack:
                current_node, path, visited = stack.pop()
                if len(path) > max_length:
                    continue
                for neighbor in self.adjacencyDict.get(current_node, []):
                    if neighbor == path[0] and len(path) >= 2:
                        if all(node >= path[0] for node in path):
                            cycle_key = tuple(path)
                            if cycle_key not in self.found_cycles:
                                self.found_cycles.add(cycle_key)
                                self.add_cycle_chain(path, is_chain=False)
                    elif neighbor not in visited:
                        visited.add(neighbor)
                        stack.append((neighbor, path + [neighbor], visited.copy()))
                        visited.remove(neighbor)

    def find_chains(self):
        max_length = self.max_chain_length
        self.found_chains = set()
        for ndd in self.NDDs:
            if max_length >= 1:
                self.add_cycle_chain([ndd.id], is_chain=True)

            for neighbor in self.adjacencyDict.get(ndd.id, []):
                stack = [(neighbor, [ndd.id, neighbor], set([neighbor]))]
                while stack:
                    current_node, path, visited = stack.pop()
                    chain_length = len(path)
                    if chain_length > max_length:
                        continue
                    self.add_cycle_chain(path, is_chain=True)
                    if chain_length < max_length:
                        for next_neighbor in self.adjacencyDict.get(current_node, []):
                            if next_neighbor not in visited:
                                visited.add(next_neighbor)
                                stack.append(
                                    (
                                        next_neighbor,
                                        path + [next_neighbor],
                                        visited.copy(),
                                    )
                                )
                                visited.remove(next_neighbor)

    def add_cycle_chain(self, nodes, is_chain):
        total_score = 0
        nbBA = 0
        num_nodes = len(nodes)
        if is_chain:
            num_arcs = num_nodes - 1
        else:
            num_arcs = num_nodes

        for i in range(num_arcs):
            from_node = nodes[i]
            to_node = nodes[(i + 1) % num_nodes]
            total_score += self.scoresDict.get((from_node, to_node), 0)

            if is_chain:
                for previous_node in nodes[:i]:
                    if self.compatibilityMatrix[from_node][previous_node] == 1:
                        nbBA += 1
            else:
                for previous_node in nodes[:i]:
                    if (
                        previous_node != nodes[i - 1]
                        and self.compatibilityMatrix[from_node][previous_node] == 1
                    ):
                        nbBA += 1
        cycle_chain = CycleChain(
            id=len(self.cyclechains),
            size=num_nodes,
            idX=nodes.copy(),
            nbBA=nbBA,
            score=total_score,
            isChain=1 if is_chain else 0,
        )
        self.cyclechains.append(cycle_chain)

    def generate_objectives(self):
        # Objective 1: Maximize the total number of transplants
        self.objectives.append(
            {
                "name": "Maximize Total Transplants",
                "sense": GRB.MAXIMIZE,
            }
        )

        # Objectives for cycles and chains sizes
        max_length = max(self.max_cycle_length, self.max_chain_length)
        for length in range(max_length, 2, -1):
            self.objectives.append(
                {
                    "name": f"Minimize Number of Cycles and Chains of Size {length}",
                    "sense": GRB.MINIMIZE,
                    "size": length,
                }
            )

        # Objective for back arcs
        self.objectives.append(
            {
                "name": "Maximize Number of Back Arcs",
                "sense": GRB.MAXIMIZE,
            }
        )

        # Objective for total score
        self.objectives.append(
            {
                "name": "Maximize Total Score/Weight",
                "sense": GRB.MAXIMIZE,
            }
        )

    def printCyclesChains(self):
        num_cycles = 0
        num_chains = 0
        for cyclechain in self.cyclechains:
            if cyclechain.isChain == 1:
                num_chains += 1
            else:
                num_cycles += 1
        print(f"Number of Cycles: {num_cycles}")
        print(f"Number of Chains: {num_chains}\n")

        for cyclechain in self.cyclechains:
            cyclechain_type = "Chain" if cyclechain.isChain == 1 else "Cycle"

            print(
                f"{cyclechain_type:<7} {cyclechain.id:<5} size {cyclechain.size:<3} ",
                end="",
            )
            for i in range(cyclechain.size - 1):
                print(f"{cyclechain.idX[i]:<4} ", end="")
            print(
                f"{cyclechain.idX[-1]:<4}  nbBackArcs {cyclechain.nbBA:<3}  score {cyclechain.score}"
            )

    def printAndWriteInfo(self, selected_cycles_chains, output_file_path):
        info = {
            "Optimal Solution Found": self.info.opt,
            "Max Cycle Length": self.max_cycle_length,
            "Max Chain Length": self.max_chain_length,
            "Total Time (s)": self.info.timeCPU[0],
            "Initialization Time (s)": self.info.timeCPU[1],
            "Total Optimization Time (s)": self.info.timeCPU[2],
        }

        # add times and failures for each objective
        for idx, obj in enumerate(self.objectives):
            obj_value = (
                self.objectiveValues[idx] if idx < len(self.objectiveValues) else None
            )
            info[f"Objective {idx + 1} ({obj['name']})"] = obj_value

            fail_count = self.fails[idx] if idx < len(self.fails) else 0
            info[f"Objective {idx + 1} Failures"] = fail_count

        info["Number of Variables"] = self.info.nbVar
        info["Number of Constraints"] = self.info.nbCons
        info["Number of Non-Zeros"] = self.info.nbNZ

        with open(output_file_path, "w") as f:
            if not selected_cycles_chains:
                f.write("No feasible solution found.\n")
                print("No feasible solution found.")
            else:
                f.write("Solution:\n")
                for idx, cc in enumerate(selected_cycles_chains):
                    f.write(f"{idx + 1}:\n")
                    f.write(f"Type: {'Chain' if cc.isChain else 'Cycle'}\n")
                    f.write(f"Size: {len(cc.idX)}\n")
                    f.write(f"Nodes: {', '.join(map(str, cc.idX))}\n")
                    f.write(f"Number of Back Arcs: {cc.nbBA}\n")
                    f.write(f"Score: {cc.score}\n")
                    f.write("-------------------------------------------\n")
                total_score = sum(cc.score for cc in selected_cycles_chains)
                f.write(f"Total score: {total_score}\n")
                print(f"Total score: {total_score}")

            f.write("\nOptimization information:\n")

            for label, value in info.items():
                f.write(f"{label}: {value}\n")
                print(f"{label}: {value}")
