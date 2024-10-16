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
        self.timeCPU = [0.0] * 7  # cpu times for different phases in optimization
        self.LB = 0  # lower bound
        self.UB = 0  # upper bound
        self.contUB = 0.0  # continuous upper bound
        self.nbCons = 0  # number of constraints
        self.nbVar = 0  # number of variables
        self.nbNZ = 0  # number of non-zeros
        self.contUB2 = 0.0  # second continuous upper bound
        self.nbCons2 = 0  # second number of constraints
        self.nbVar2 = 0  # second number of variables
        self.nbNZ2 = 0  # second number of non-zeros


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
    def __init__(self):
        self.NDDs = []
        self.pairs = []
        self.cyclechains = []
        self.maxId = 0
        self.nbPairs = 0
        self.nbNDDs = 0
        self.idToIdxA = {}
        self.idToIdxP = {}
        self.comp = []
        self.scores = []
        self.info = Info()
        self.isActivated = []
        self.objs = []
        self.RC = []
        self.tObjs = []
        self.fails = [0] * 4
        self.types = [0] * 6

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

        self.comp = [
            [-1 for _ in range(self.maxId + 1)] for _ in range(self.maxId + 1)
        ]  # initialize
        self.scores = [
            [-1 for _ in range(self.maxId + 1)] for _ in range(self.maxId + 1)
        ]  # initialize

        # fill compatibility and score matrices based on arcs
        for arc in data["arcs"]:
            donor_id = arc["donor_id"]
            patient_id = arc["patient_id"]
            weight = arc["weight"]

            self.comp[donor_id][patient_id] = 1
            self.scores[donor_id][patient_id] = weight

        # link targets and scores to pairs and NDDs
        for pair in self.pairs:
            pair.nbEdges = sum(
                [1 for arc in data["arcs"] if arc["donor_id"] == pair.source]
            )
            pair.targets = [
                arc["patient_id"]
                for arc in data["arcs"]
                if arc["donor_id"] == pair.source
            ]
            pair.scores = [
                arc["weight"] for arc in data["arcs"] if arc["donor_id"] == pair.source
            ]

        for ndd in self.NDDs:
            ndd.nbEdges = sum([1 for arc in data["arcs"] if arc["donor_id"] == ndd.id])
            ndd.targets = [
                arc["patient_id"] for arc in data["arcs"] if arc["donor_id"] == ndd.id
            ]
            ndd.scores = [
                arc["weight"] for arc in data["arcs"] if arc["donor_id"] == ndd.id
            ]

        # create cycles
        for i in range(self.maxId + 1):
            if i in self.idToIdxP:
                hbP = [False] * (self.maxId + 1)
                for idx in self.idToIdxP[i]:
                    for k, tId in enumerate(self.pairs[idx].targets):
                        if not hbP[tId]:
                            if self.comp[tId][i] == 1 and i <= tId:
                                cycle = CycleChain(
                                    id=len(self.cyclechains),
                                    size=2,
                                    idX=[i, tId],
                                    nbBA=0,
                                    score=self.scores[i][tId] + self.scores[tId][i],
                                )
                                self.cyclechains.append(cycle)
                                self.types[0] += 1

                            hbP[tId] = True
                            hbP2 = [False] * (self.maxId + 1)
                            for idx2 in self.idToIdxP.get(tId, []):
                                for m, tId2 in enumerate(self.pairs[idx2].targets):
                                    if (
                                        not hbP2[tId2]
                                        and self.comp[tId2][i] == 1
                                        and i <= tId
                                        and i <= tId2
                                    ):
                                        cycle = CycleChain(
                                            id=len(self.cyclechains),
                                            size=3,
                                            idX=[i, tId, tId2],
                                            nbBA=0,
                                            score=self.scores[i][tId]
                                            + self.scores[tId][tId2]
                                            + self.scores[tId2][i],
                                        )
                                        if self.comp[i][tId2] == 1:
                                            cycle.nbBA += 1
                                        if self.comp[tId][i] == 1:
                                            cycle.nbBA += 1
                                        if self.comp[tId2][tId] == 1:
                                            cycle.nbBA += 1
                                        self.cyclechains.append(cycle)
                                        self.types[1] += 1
                                    hbP2[tId2] = True
                                    # TODO: extend to more than 3 cycles

        # create chains (NDDs initiating chains)
        for ndd in self.NDDs:
            chain1 = CycleChain(
                id=len(self.cyclechains),
                size=1,
                idX=[ndd.id],
                nbBA=0,
                score=0,
                isChain=1,
            )
            self.cyclechains.append(chain1)
            self.types[2] += 1

            for j, tId in enumerate(ndd.targets):
                chain2 = CycleChain(
                    id=len(self.cyclechains),
                    size=2,
                    idX=[ndd.id, tId],
                    nbBA=0,
                    score=self.scores[ndd.id][tId],
                    isChain=1,
                )
                self.cyclechains.append(chain2)
                self.types[3] += 1
                hbP2 = [False] * (self.maxId + 1)
                for idx2 in self.idToIdxP.get(tId, []):
                    for m, tId2 in enumerate(self.pairs[idx2].targets):
                        if not hbP2[tId2]:
                            chain3 = CycleChain(
                                id=len(self.cyclechains),
                                size=3,
                                idX=[ndd.id, tId, tId2],
                                nbBA=0,
                                score=self.scores[ndd.id][tId] + self.scores[tId][tId2],
                                isChain=1,
                            )
                            if self.comp[ndd.id][tId2] == 1:
                                chain3.nbBA += 1
                            if self.comp[tId2][tId] == 1:
                                chain3.nbBA += 1
                            self.cyclechains.append(chain3)
                            self.types[4] += 1
                            hbP3 = [False] * (self.maxId + 1)
                            for idx3 in self.idToIdxP.get(tId2, []):
                                for o, tId3 in enumerate(self.pairs[idx3].targets):
                                    if not hbP3[tId3] and tId != tId3:
                                        chain4 = CycleChain(
                                            id=len(self.cyclechains),
                                            size=4,
                                            idX=[ndd.id, tId, tId2, tId3],
                                            nbBA=0,
                                            score=self.scores[ndd.id][tId]
                                            + self.scores[tId][tId2]
                                            + self.scores[tId2][tId3],
                                            isChain=1,
                                        )
                                        if self.comp[ndd.id][tId2] == 1:
                                            chain4.nbBA += 1
                                        if self.comp[ndd.id][tId3] == 1:
                                            chain4.nbBA += 1
                                        if self.comp[tId2][tId] == 1:
                                            chain4.nbBA += 1
                                        if self.comp[tId][tId3] == 1:
                                            chain4.nbBA += 1
                                        if self.comp[tId3][tId] == 1:
                                            chain4.nbBA += 1
                                        self.cyclechains.append(chain4)
                                        self.types[5] += 1
                                    hbP3[tId3] = True
                        hbP2[tId2] = True

    def printProb(self):
        print(f"Instance")
        for ndd in self.NDDs:
            ndd.print()
        for pair in self.pairs:
            pair.print()
        print(f"Total number of Cycles and Chains: {len(self.cyclechains)}")

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

    def printAndWriteInfo(allo, selected_cycles_chains, output_file_path):
        info = {
            "Optimal Solution Found": allo.info.opt,
            "Total Time (s)": allo.info.timeCPU[0],
            "Initialization Time (s)": allo.info.timeCPU[1],
            "CycleLP Time (s)": allo.info.timeCPU[2],
            "Cycle Deactivation Time (s)": allo.info.timeCPU[3],
            "Other Time (s)": allo.info.timeCPU[4],
            "Final Optimization Time (s)": allo.info.timeCPU[5],
            "Post-Processing Time (s)": allo.info.timeCPU[6],
            "Objective 1 (Max. Cycles and Chains)": (
                allo.objs[0] if len(allo.objs) > 0 else None
            ),
            "Objective 2 (Min. Cycles and Chains of Size 4)": (
                allo.objs[1] if len(allo.objs) > 1 else None
            ),
            "Objective 3 (Min. Cycles Chains of Size 3)": (
                allo.objs[2] if len(allo.objs) > 2 else None
            ),
            "Objective 4 (Max. Number of Backarcs)": (
                allo.objs[3] if len(allo.objs) > 3 else None
            ),
            "Objective 5 (Max. Total Score/Weight)": (
                allo.objs[4] if len(allo.objs) > 4 else None
            ),
            "Number of Variables": allo.info.nbVar,
            "Number of Constraints": allo.info.nbCons,
            "Number of Non-Zeros": allo.info.nbNZ,
            "CycleLP Failures": allo.fails[0] if len(allo.fails) > 0 else None,
            "Cycle Deactivation Failures": (
                allo.fails[1] if len(allo.fails) > 1 else None
            ),
            "Post-Processing Failures": allo.fails[2] if len(allo.fails) > 2 else None,
        }

        with open(output_file_path, "w") as f:
            f.write("Solution:\n")
            for idx, cc in enumerate(selected_cycles_chains):
                f.write(f"{idx + 1}:\n")
                if cc.isChain:
                    f.write("Type: Chain\n")
                else:
                    f.write("Type: Cycle\n")
                f.write(f"Size: {len(cc.idX)}\n")
                f.write(f"Nodes: {', '.join(map(str, cc.idX))}\n")
                f.write(f"Number of Back Arcs: {cc.nbBA}\n")
                f.write(f"Score: {cc.score}\n")
                f.write("-------------------------------------------\n")
            f.write(f"Total score: {sum(cc.score for cc in selected_cycles_chains)}\n")
            f.write("\nOptimization information:\n")

            for label, value in info.items():
                f.write(f"{label}: {value}\n")
                print(f"{label}: {value}")
