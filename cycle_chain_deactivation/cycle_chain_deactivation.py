import time
import math
import gurobipy as gp
from gurobipy import GRB
from allocation import Allocation

EPSILON = 0.001
TIMEOUT = 300

def run_cycle_chain_deactivation(data):
    init_time_model_cpu = time.time()

    allo = Allocation()
    allo.load(data)
    # allo.printProb()
    allo.info.timeCPU[1] = time.time() - init_time_model_cpu

    allo.isActivated = [1] * len(allo.cyclechains)

    for i in range(4):
        allo.fails.append(0)
        cycleLP(allo, init_time_model_cpu)
        if i == 0 or i == 3:
            sol = math.floor(allo.tObjs[i] + EPSILON)
        else:
            sol = math.ceil(allo.tObjs[i] - EPSILON)
        

        # deactivation loop
        while True:
            if time.time() - init_time_model_cpu > TIMEOUT:
                allo.objs[i] = -1
                allo.info['opt'] = False
                break

            print(f"Sol {i} is at {sol}")

            # deactivation logic
            for j in range(len(allo.cyclechains)):
                if i == 0 or i == 3:
                    if allo.isActivated[j] >= 0:
                        if allo.tObjs[i] + allo.RC[j] + EPSILON < sol:
                            allo.isActivated[j] = 0
                        else:
                            allo.isActivated[j] = 1
                else:
                    if allo.isActivated[j] >= 0:
                        if allo.tObjs[i] + allo.RC[j] - EPSILON > sol:
                            allo.isActivated[j] = 0
                        else:
                            allo.isActivated[j] = 1

            # check if cycle conditions are met
            if cycleILP(allo, init_time_model_cpu) == 0 and allo.objs[i] == sol:
                for j in range(len(allo.cyclechains)):
                    if allo.isActivated[j] == 0:
                        allo.isActivated[j] = -1
                break
            else:
                allo.objs.append(sol)
                if time.time() - init_time_model_cpu < TIMEOUT:
                    allo.fails[i] += 1
                    if i == 0 or i == 3:
                        sol -= 1
                    else:
                        sol += 1

        # record CPU time
        allo.info.timeCPU[i + 2] = time.time() - init_time_model_cpu
        for j in range(i + 1):
            allo.info.timeCPU[i + 2] -= allo.info.timeCPU[j + 1]

    cycleILP(allo, init_time_model_cpu)
    allo.info.timeCPU[6] = time.time() - init_time_model_cpu
    for j in range(5):
        allo.info.timeCPU[6] -= allo.info.timeCPU[j + 1]

    allo.info.timeCPU[0] = time.time() - init_time_model_cpu

    allo.printInfo()


def cycleLP(allo, init_time_model_cpu):
    try:
        model = gp.Model("cycleLP")

        objFun1 = gp.LinExpr(0)
        objFun2 = gp.LinExpr(0)
        objFun3 = gp.LinExpr(0)
        objFun4 = gp.LinExpr(0)
        objFun5 = gp.LinExpr(0)

        # ===== something is going wrong here ====== -->
        isCycleUsed = [None] * len(allo.cyclechains)
        isPatientUsed = [gp.LinExpr(0)] * (allo.maxId + 1)
        isPatientIdUsed = [False] * (allo.maxId + 1)

        # Initialization
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                isCycleUsed[i] = model.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS)

        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                    for j in allo.cyclechains[i].idX:
                        if j not in allo.idToIdxA:
                            isPatientUsed[j] += isCycleUsed[i]
                            isPatientIdUsed[j] = True
                    objFun1 += len(allo.cyclechains[i].idX) * isCycleUsed[i]
                    if len(allo.cyclechains[i].idX) == 4:
                        objFun2 += isCycleUsed[i]
                    if len(allo.cyclechains[i].idX) == 3:
                        objFun3 += isCycleUsed[i]
                    objFun4 += allo.cyclechains[i].nbBA * isCycleUsed[i]
                    objFun5 += allo.cyclechains[i].score * isCycleUsed[i]
        # <-- ===== something is going wrong here ======
                

        # each patient can only be used once
        for i in range(allo.maxId + 1):
            if isPatientIdUsed[i]:
                model.addConstr(isPatientUsed[i] <= 1)

        # setting objective based on number of objectives
        if len(allo.objs) == 0:
            model.setObjective(objFun1, GRB.MAXIMIZE)
        elif len(allo.objs) == 1:
            model.addConstr(objFun1 == allo.objs[0])
            model.setObjective(objFun2, GRB.MINIMIZE)
        elif len(allo.objs) == 2:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.setObjective(objFun3, GRB.MINIMIZE)
        elif len(allo.objs) == 3:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            model.setObjective(objFun4, GRB.MAXIMIZE)
        elif len(allo.objs) == 4:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            model.addConstr(objFun4 == allo.objs[3])
            model.setObjective(objFun5, GRB.MAXIMIZE)

        # setting Gurobi parameters
        model.setParam('TimeLimit', TIMEOUT - (time.time() - init_time_model_cpu))
        model.setParam('Threads', 1)
        model.setParam('Method', 2)
        model.setParam('Crossover', 0)
        model.setParam('MIPGap', 0)
        model.optimize()

        # if solution found
        allo.tObjs.append(model.objVal)

        # fill solution
        allo.RC = [0] * len(allo.cyclechains)
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                if isCycleUsed[i].X < EPSILON:
                    allo.RC[i] = isCycleUsed[i].RC
                else:
                    allo.RC[i] = 0.0
                if i < 0:
                    print(f"{isCycleUsed[i].X} ({isCycleUsed[i].RC}) x ")
                    allo.cyclechains[i].print()
        
    except gp.GurobiError as e:
        print(f"Error code = {e.errno}")
        print(e.message)
    except Exception as e:
        print("Exception during optimization:", e)


def cycleILP(allo, init_time_model_cpu):
    try:
        model = gp.Model("cycle")

        objFun1 = gp.LinExpr(0)
        objFun2 = gp.LinExpr(0)
        objFun3 = gp.LinExpr(0)
        objFun4 = gp.LinExpr(0)
        objFun5 = gp.LinExpr(0)

        # ===== something is going wrong here ====== -->
        isCycleUsed = [None] * len(allo.cyclechains)
        isPatientUsed = [gp.LinExpr(0)] * (allo.maxId + 1)
        isPatientIdUsed = [False] * (allo.maxId + 1)

        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                isCycleUsed[i] = model.addVar(lb=0, ub=1, vtype=GRB.INTEGER)

        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                for j in allo.cyclechains[i].idX:
                    if j not in allo.idToIdxA:
                        isPatientUsed[j] += isCycleUsed[i]
                        isPatientIdUsed[j] = True
                objFun1 += len(allo.cyclechains[i].idX) * isCycleUsed[i]
                if len(allo.cyclechains[i].idX) == 4:
                    objFun2 += isCycleUsed[i]
                if len(allo.cyclechains[i].idX) == 3:
                    objFun3 += isCycleUsed[i]
                objFun4 += allo.cyclechains[i].nbBA * isCycleUsed[i]
                objFun5 += allo.cyclechains[i].score * isCycleUsed[i]
        # <-- ===== something is going wrong here ======
        
        for i in range(allo.maxId + 1):
            if isPatientIdUsed[i]:
                model.addConstr(isPatientUsed[i] <= 1)
        
        if len(allo.objs) == 0:
            model.setObjective(objFun1, GRB.MAXIMIZE)
        elif len(allo.objs) == 1:
            model.addConstr(objFun1 == allo.objs[0])
            # model.setObjective(objFun2, GRB.MINIMIZE)
        elif len(allo.objs) == 2:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            # model.setObjective(objFun3, GRB.MINIMIZE)
        elif len(allo.objs) == 3:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            # model.setObjective(objFun4, GRB.MAXIMIZE)
        elif len(allo.objs) == 4:
            model.addConstr(objFun1 == allo.objs[0])
            model.addConstr(objFun2 == allo.objs[1])
            model.addConstr(objFun3 == allo.objs[2])
            model.addConstr(objFun4 == allo.objs[3])
            model.setObjective(objFun5, GRB.MAXIMIZE)

        model.setParam('TimeLimit', TIMEOUT - (time.time() - init_time_model_cpu))
        model.setParam('Threads', 1)
        model.setParam('Method', 2)
        model.setParam('MIPGap', 0)
        model.optimize()

        allo.info.UB = math.ceil(model.ObjBound - EPSILON)
        allo.info.opt = False

        allo.info.nbVar = model.NumVars
        allo.info.nbCons = model.NumConstrs
        allo.info.nbNZ = model.NumNZs

        # if no solution found
        if model.SolCount < 1:
            print("Failed to optimize ILP.")
            allo.info.LB = 0
            return -1
        
        # if solution found
        allo.info.LB = math.ceil(model.ObjVal - EPSILON)
        if allo.info.LB == allo.info.UB:
            allo.info.opt = True
        allo.objs.append(math.ceil(model.ObjVal - EPSILON))

        # filling solution
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1:
                if isCycleUsed[i].X > EPSILON:
                    print(f"{isCycleUsed[i].X} x ")
                    allo.cyclechains[i].print()

    except gp.GurobiError as e:
        print(f"Error code = {e.errno}")
        print(e.message)
    except Exception as e:
        print("Exception during optimization:", e)

    return 0
