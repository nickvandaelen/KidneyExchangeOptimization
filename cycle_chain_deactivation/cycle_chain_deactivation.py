import time
import math
import gurobipy as gp
from gurobipy import GRB
from allocation import Allocation

EPSILON = 0.001
TIMEOUT = 300
MAX_ITERATIONS = 1000


def run_cycle_chain_deactivation(data, output_file_path):
    init_time_model_cpu = time.time()

    allo = Allocation()
    allo.load(data)
    allo.info.timeCPU[1] = time.time() - init_time_model_cpu

    allo.isActivated = [1] * len(allo.cyclechains)
    allo.objs = [None] * 5
    allo.tObjs = [0] * 5
    allo.fails = [0] * 4

    for i in range(4):
        cycleLP(allo, init_time_model_cpu, i)
        if i == 0 or i == 3:
            sol = math.floor(allo.tObjs[i] + EPSILON)
        else:
            sol = math.ceil(allo.tObjs[i] - EPSILON)

        # deactivation loop
        iteration = 0
        while True:
            if time.time() - init_time_model_cpu > TIMEOUT:
                allo.objs[i] = -1
                allo.info.opt = False
                break

            print(f"Iteration {iteration}, Sol {i} is at {sol}")

            # deactivation logic
            for j in range(len(allo.cyclechains)):
                if allo.isActivated[j] >= 0:
                    if i == 0 or i == 3:
                        if allo.tObjs[i] + allo.RC[j] + EPSILON < sol:
                            allo.isActivated[j] = 0
                        else:
                            allo.isActivated[j] = 1
                    else:
                        if allo.tObjs[i] + allo.RC[j] - EPSILON > sol:
                            allo.isActivated[j] = 0
                        else:
                            allo.isActivated[j] = 1

            # solve the ILP model for the current objective
            obj_val, _ = cycleILP(allo, init_time_model_cpu, i)
            print(
                f"Iteration {iteration}, Objective {i}, Sol = {sol}, ObjVal = {obj_val}"
            )

            if obj_val == -1:
                print(f"Failed to optimize ILP for objective {i}.")
                allo.objs[i] = -1
                break
            elif obj_val == sol:
                allo.objs[i] = obj_val
                # deactivate cycles/chains permanently
                for j in range(len(allo.cyclechains)):
                    if allo.isActivated[j] == 0:
                        allo.isActivated[j] = -1
                break
            else:
                allo.objs[i] = sol
                if time.time() - init_time_model_cpu < TIMEOUT:
                    allo.fails[i] += 1
                    # adjust sol based on the objective
                    if i == 0 or i == 3:
                        sol -= 1
                    else:
                        sol += 1
            iteration += 1
            if iteration >= MAX_ITERATIONS:
                print(f"Reached maximum iterations for objective {i}")
                break
        allo.info.timeCPU[i + 2] = time.time() - init_time_model_cpu
        for j in range(i + 1):
            allo.info.timeCPU[i + 2] -= allo.info.timeCPU[j + 1]

    obj_val, selected_cycles_chains = cycleILP(allo, init_time_model_cpu, 4)
    allo.info.timeCPU[6] = time.time() - init_time_model_cpu
    for j in range(5):
        allo.info.timeCPU[6] -= allo.info.timeCPU[j + 1]

    allo.info.timeCPU[0] = time.time() - init_time_model_cpu

    allo.printAndWriteInfo(selected_cycles_chains, output_file_path)


def cycleLP(allo, init_time_model_cpu, objective_index):
    try:
        model = gp.Model("cycleLP")

        objFun1 = gp.LinExpr(0)
        objFun2 = gp.LinExpr(0)
        objFun3 = gp.LinExpr(0)
        objFun4 = gp.LinExpr(0)
        objFun5 = gp.LinExpr(0)

        isCycleUsed = [None] * len(allo.cyclechains)
        isPatientUsed = [gp.LinExpr(0) for _ in range(allo.maxId + 1)]
        isPatientIdUsed = [False] * (allo.maxId + 1)

        # add variables for activated cycles/chains
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                isCycleUsed[i] = model.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS)

        # build objective functions and constraints
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                for j in allo.cyclechains[i].idX:
                    isPatientUsed[j] += isCycleUsed[i]
                    isPatientIdUsed[j] = True
                if allo.cyclechains[i].isChain:
                    objFun1 += (len(allo.cyclechains[i].idX) - 1) * isCycleUsed[i]
                else:
                    objFun1 += len(allo.cyclechains[i].idX) * isCycleUsed[i]
                if len(allo.cyclechains[i].idX) == 4:
                    objFun2 += isCycleUsed[i]
                if len(allo.cyclechains[i].idX) == 3:
                    objFun3 += isCycleUsed[i]
                objFun4 += allo.cyclechains[i].nbBA * isCycleUsed[i]
                objFun5 += allo.cyclechains[i].score * isCycleUsed[i]

        # each patient can be used at most once
        for i in range(allo.maxId + 1):
            if isPatientIdUsed[i]:
                model.addConstr(isPatientUsed[i] <= 1)

        # setting objective based on current objective index
        if objective_index == 0:
            model.setObjective(objFun1, GRB.MAXIMIZE)
        elif objective_index == 1:
            model.addConstr(objFun1 == allo.objs[0])
            model.setObjective(objFun2, GRB.MINIMIZE)
        elif objective_index == 2:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.setObjective(objFun3, GRB.MINIMIZE)
        elif objective_index == 3:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            model.setObjective(objFun4, GRB.MAXIMIZE)
        elif objective_index == 4:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            model.addConstr(objFun4 == allo.objs[3])
            model.setObjective(objFun5, GRB.MAXIMIZE)

        # setting Gurobi parameters
        model.setParam("TimeLimit", TIMEOUT - (time.time() - init_time_model_cpu))
        model.setParam("Method", -1)
        model.setParam("Crossover", 0)
        model.setParam("MIPGap", 0)
        model.optimize()

        # store objective value
        allo.tObjs[objective_index] = model.ObjVal

        # compute reduced costs
        allo.RC = [0] * len(allo.cyclechains)
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                if isCycleUsed[i].X < EPSILON:
                    allo.RC[i] = isCycleUsed[i].RC
                else:
                    allo.RC[i] = 0.0

    except gp.GurobiError as e:
        print(f"Error code = {e.errno}")
        print(e.message)
    except Exception as e:
        print("Exception during optimization:", e)


def cycleILP(allo, init_time_model_cpu, objective_index):
    try:
        model = gp.Model("cycle")

        objFun1 = gp.LinExpr(0)
        objFun2 = gp.LinExpr(0)
        objFun3 = gp.LinExpr(0)
        objFun4 = gp.LinExpr(0)
        objFun5 = gp.LinExpr(0)

        isCycleUsed = [None] * len(allo.cyclechains)
        isPatientUsed = [gp.LinExpr(0) for _ in range(allo.maxId + 1)]
        isPatientIdUsed = [False] * (allo.maxId + 1)

        # add variables for activated cycles/chains
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                isCycleUsed[i] = model.addVar(lb=0, ub=1, vtype=GRB.BINARY)

        # build objective functions and constraints
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                for j in allo.cyclechains[i].idX:
                    isPatientUsed[j] += isCycleUsed[i]
                    isPatientIdUsed[j] = True
                if allo.cyclechains[i].isChain:
                    objFun1 += (len(allo.cyclechains[i].idX) - 1) * isCycleUsed[i]
                else:
                    objFun1 += len(allo.cyclechains[i].idX) * isCycleUsed[i]
                if len(allo.cyclechains[i].idX) == 4:
                    objFun2 += isCycleUsed[i]
                if len(allo.cyclechains[i].idX) == 3:
                    objFun3 += isCycleUsed[i]
                objFun4 += allo.cyclechains[i].nbBA * isCycleUsed[i]
                objFun5 += allo.cyclechains[i].score * isCycleUsed[i]

        # each patient can be used at most once
        for i in range(allo.maxId + 1):
            if isPatientIdUsed[i]:
                model.addConstr(isPatientUsed[i] <= 1)

        # set the objective and constraints based on the current objective index
        if objective_index == 0:
            model.setObjective(objFun1, GRB.MAXIMIZE)
        elif objective_index == 1:
            model.addConstr(objFun1 == allo.objs[0])
            model.setObjective(objFun2, GRB.MINIMIZE)
        elif objective_index == 2:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.setObjective(objFun3, GRB.MINIMIZE)
        elif objective_index == 3:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            model.setObjective(objFun4, GRB.MAXIMIZE)
        elif objective_index == 4:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            model.addConstr(objFun4 == allo.objs[3])
            model.setObjective(objFun5, GRB.MAXIMIZE)

        # Set Gurobi parameters
        model.setParam("TimeLimit", TIMEOUT - (time.time() - init_time_model_cpu))
        model.setParam("Method", -1)
        model.setParam("MIPGap", 0)
        model.optimize()

        # update allocation information
        allo.info.UB = math.ceil(model.ObjBound - EPSILON)
        allo.info.opt = False

        allo.info.nbVar = model.NumVars
        allo.info.nbCons = model.NumConstrs
        allo.info.nbNZ = model.NumNZs

        # check if no solution found
        if model.SolCount < 1:
            print("Failed to optimize ILP.")
            allo.info.LB = 0
            return -1

        # store the objective value
        obj_val = math.ceil(model.ObjVal - EPSILON)
        allo.objs[objective_index] = obj_val

        # update optimality status
        allo.info.LB = obj_val
        if allo.info.LB == allo.info.UB:
            allo.info.opt = True

        # collect and print the solution
        selected_cycles_chains = []
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1 and isCycleUsed[i].X > EPSILON:
                selected_cycles_chains.append(allo.cyclechains[i])
                print(f"{isCycleUsed[i].X} x ")
                allo.cyclechains[i].print()

        return obj_val, selected_cycles_chains

    except gp.GurobiError as e:
        print(f"Error code = {e.errno}")
        print(e.message)
        return -1, []
    except Exception as e:
        print("Exception during optimization:", e)
        return -1, []
