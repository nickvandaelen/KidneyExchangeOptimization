import time
import math
import gurobipy as gp
from gurobipy import GRB
from allocation_generalized import Allocation

EPSILON = 0.001
TIMEOUT = 7200
MAX_ITERATIONS = 2000


def run_cycle_chain_deactivation(
    data, output_file_path, max_cycle_length=3, max_chain_length=4
):
    start_time = time.time()

    allo = Allocation(max_cycle_length, max_chain_length)
    allo.load(data)
    initialization_time = time.time() - start_time

    start_optimization_time = time.time()

    num_objectives = len(allo.objectives)
    allo.isActivated = [1] * len(allo.cyclechains)
    allo.objectiveValues = [None] * num_objectives
    allo.temporaryObjectiveValues = [0] * num_objectives
    allo.fails = [0] * num_objectives

    for i in range(num_objectives - 1):
        cycleLP(allo, start_optimization_time, i)
        current_obj = allo.objectives[i]
        if current_obj["sense"] == GRB.MAXIMIZE:
            sol = math.floor(allo.temporaryObjectiveValues[i] + EPSILON)
        else:
            sol = math.ceil(allo.temporaryObjectiveValues[i] - EPSILON)

        # deactivation loop
        iteration = 0
        while True:
            if time.time() - start_optimization_time > TIMEOUT:
                allo.objectiveValues[i] = -1
                allo.info.opt = False
                break

            print(
                f"Iteration {iteration}, Objective {i} ({current_obj['name']}) is at {sol}"
            )

            # deactivation logic
            for j in range(len(allo.cyclechains)):
                if allo.isActivated[j] >= 0:
                    if current_obj["sense"] == GRB.MAXIMIZE:
                        if (
                            allo.temporaryObjectiveValues[i] + allo.RC[j] + EPSILON
                            < sol
                        ):
                            allo.isActivated[j] = 0
                        else:
                            allo.isActivated[j] = 1
                    else:
                        if (
                            allo.temporaryObjectiveValues[i] + allo.RC[j] - EPSILON
                            > sol
                        ):
                            allo.isActivated[j] = 0
                        else:
                            allo.isActivated[j] = 1

            # solve the ILP model for the current objective
            obj_val, _ = cycleILP(allo, start_optimization_time, i)
            print(
                f"Iteration {iteration}, Objective {i}, Sol = {sol}, ObjVal = {obj_val}"
            )

            if obj_val == -1:
                print(f"Failed to optimize ILP for objective {i}.")
                allo.objectiveValues[i] = -1
                break
            elif obj_val == sol:
                allo.objectiveValues[i] = obj_val
                # deactivate cycles/chains permanently
                for j in range(len(allo.cyclechains)):
                    if allo.isActivated[j] == 0:
                        allo.isActivated[j] = -1
                break
            else:
                allo.objectiveValues[i] = sol
                if time.time() - start_optimization_time < TIMEOUT:
                    allo.fails[i] += 1
                    # adjust sol based on the objective
                    if current_obj["sense"] == GRB.MAXIMIZE:
                        sol -= 1
                    else:
                        sol += 1
            iteration += 1
            if iteration >= MAX_ITERATIONS:
                print(f"Reached maximum iterations for objective {i}")
                break

    obj_val, selected_cycles_chains = cycleILP(
        allo, start_optimization_time, num_objectives - 1
    )
    if obj_val == -1:
        print("No feasible solution found in the final optimization.")
        allo.info.opt = False
        selected_cycles_chains = []

    end_time = time.time()
    total_optimization_time = end_time - start_optimization_time
    total_time = end_time - start_time

    # store relevant times
    allo.info.timeCPU.append(total_time)
    allo.info.timeCPU.append(initialization_time)
    allo.info.timeCPU.append(total_optimization_time)

    allo.printAndWriteInfo(selected_cycles_chains, output_file_path)


def cycleLP(allo, start_optimization_time, objective_index):
    try:
        model = gp.Model("cycleLP")

        isCycleUsed = [None] * len(allo.cyclechains)
        isPatientUsed = [gp.LinExpr(0) for _ in range(allo.maxId + 1)]
        isPatientIdUsed = [False] * (allo.maxId + 1)

        # initialize objective functions
        objFunTransplants = gp.LinExpr(0)
        objFunBackArcs = gp.LinExpr(0)
        objFunScore = gp.LinExpr(0)
        objFunSizeDict = {}

        max_length = max(allo.max_cycle_length, allo.max_chain_length)
        for length in range(max_length, 1, -1):
            objFunSizeDict[length] = gp.LinExpr(0)

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

                size = len(allo.cyclechains[i].idX)

                if allo.cyclechains[i].isChain:
                    objFunTransplants += (size - 1) * isCycleUsed[i]
                else:
                    objFunTransplants += size * isCycleUsed[i]

                if size in objFunSizeDict:
                    objFunSizeDict[size] += isCycleUsed[i]

                objFunBackArcs += allo.cyclechains[i].nbBA * isCycleUsed[i]
                objFunScore += allo.cyclechains[i].score * isCycleUsed[i]

        # each patient can be used at most once
        for i in range(allo.maxId + 1):
            if isPatientIdUsed[i]:
                model.addConstr(isPatientUsed[i] <= 1)

        # add constraints for previous objectives
        for idx in range(objective_index):
            prev_obj = allo.objectives[idx]
            if idx == 0:
                model.addConstr(objFunTransplants == allo.objectiveValues[idx])
            elif "size" in prev_obj:
                model.addConstr(
                    objFunSizeDict[prev_obj["size"]] == allo.objectiveValues[idx]
                )
            elif prev_obj["name"] == "Maximize Number of Back Arcs":
                model.addConstr(objFunBackArcs == allo.objectiveValues[idx])
            elif prev_obj["name"] == "Maximize Total Score/Weight":
                model.addConstr(objFunScore == allo.objectiveValues[idx])

        # set the current objective
        current_objective = allo.objectives[objective_index]
        if current_objective["name"] == "Maximize Total Transplants":
            model.setObjective(objFunTransplants, GRB.MAXIMIZE)
        elif "size" in current_objective:
            model.setObjective(
                objFunSizeDict[current_objective["size"]], current_objective["sense"]
            )
        elif current_objective["name"] == "Maximize Number of Back Arcs":
            model.setObjective(objFunBackArcs, GRB.MAXIMIZE)
        elif current_objective["name"] == "Maximize Total Score/Weight":
            model.setObjective(objFunScore, GRB.MAXIMIZE)

        # setting Gurobi parameters
        model.setParam("TimeLimit", TIMEOUT - (time.time() - start_optimization_time))
        model.setParam("MIPGap", 0)
        model.optimize()

        # store objective value
        allo.temporaryObjectiveValues[objective_index] = model.ObjVal

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


def cycleILP(allo, start_optimization_time, objective_index):
    try:
        model = gp.Model("cycleILP")

        isCycleUsed = [None] * len(allo.cyclechains)
        isPatientUsed = [gp.LinExpr(0) for _ in range(allo.maxId + 1)]
        isPatientIdUsed = [False] * (allo.maxId + 1)

        # initialize objective functions
        objFunTransplants = gp.LinExpr(0)
        objFunBackArcs = gp.LinExpr(0)
        objFunScore = gp.LinExpr(0)
        objFunSizeDict = {}

        max_length = max(allo.max_cycle_length, allo.max_chain_length)
        for length in range(max_length, 1, -1):
            objFunSizeDict[length] = gp.LinExpr(0)

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

                size = len(allo.cyclechains[i].idX)

                if allo.cyclechains[i].isChain:
                    objFunTransplants += (size - 1) * isCycleUsed[i]
                else:
                    objFunTransplants += size * isCycleUsed[i]

                if size in objFunSizeDict:
                    objFunSizeDict[size] += isCycleUsed[i]

                objFunBackArcs += allo.cyclechains[i].nbBA * isCycleUsed[i]
                objFunScore += allo.cyclechains[i].score * isCycleUsed[i]

        # each patient can be used at most once
        for i in range(allo.maxId + 1):
            if isPatientIdUsed[i]:
                model.addConstr(isPatientUsed[i] <= 1)

        # add constraints for previous objectives
        total_transplants_constr = None
        for idx in range(objective_index):
            prev_obj = allo.objectives[idx]
            if idx == 0:
                total_transplants_constr = model.addConstr(
                    objFunTransplants == allo.objectiveValues[idx]
                )
            elif "size" in prev_obj:
                model.addConstr(
                    objFunSizeDict[prev_obj["size"]] == allo.objectiveValues[idx]
                )
            elif prev_obj["name"] == "Maximize Number of Back Arcs":
                model.addConstr(objFunBackArcs == allo.objectiveValues[idx])
            elif prev_obj["name"] == "Maximize Total Score/Weight":
                model.addConstr(objFunScore == allo.objectiveValues[idx])

        # set the current objective
        current_objective = allo.objectives[objective_index]
        if current_objective["name"] == "Maximize Total Transplants":
            model.setObjective(objFunTransplants, GRB.MAXIMIZE)
        elif "size" in current_objective:
            model.setObjective(
                objFunSizeDict[current_objective["size"]],
                current_objective["sense"],
            )
        elif current_objective["name"] == "Maximize Number of Back Arcs":
            model.setObjective(objFunBackArcs, GRB.MAXIMIZE)
        elif current_objective["name"] == "Maximize Total Score/Weight":
            model.setObjective(objFunScore, GRB.MAXIMIZE)

        # set Gurobi parameters
        model.setParam("TimeLimit", TIMEOUT - (time.time() - start_optimization_time))
        model.setParam("MIPGap", 0)
        model.optimize()

        # check if optimal number of transplants is causing infeasible solution
        if model.Status == GRB.INFEASIBLE:
            print(
                "Model is infeasible. Checking if total transplants constraint is the cause."
            )
            model.remove(total_transplants_constr)
            model.update()
            model.optimize()

            if model.Status == GRB.INFEASIBLE:
                print("Model is still infeasible without total transplants constraint.")
                allo.objectiveValues[objective_index] = -1
                allo.info.opt = False
                return -1, []
            else:
                print("Infeasibility was due to total transplants constraint.")
                # Proceed with iterative reduction
                optimalNumTransplants = allo.objectiveValues[0]
                iteration = 1

            while True:
                required_transplants = optimalNumTransplants - iteration
                if required_transplants < 0:
                    print(
                        "No feasible solution found within acceptable number of transplants."
                    )
                    allo.objectiveValues[objective_index] = -1
                    allo.info.opt = False
                    return -1, []

                total_transplants_constr = model.addConstr(
                    objFunTransplants == required_transplants
                )
                model.update()
                model.optimize()

                if model.Status != GRB.INFEASIBLE:
                    allo.objectiveValues[0] = required_transplants
                    break
                else:
                    model.remove(total_transplants_constr)
                    model.update()
                    iteration += 1

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
            return -1, []

        # store the objective value
        obj_val = math.ceil(model.ObjVal - EPSILON)
        allo.objectiveValues[objective_index] = obj_val

        # update optimality status
        allo.info.LB = obj_val
        if allo.info.LB == allo.info.UB:
            allo.info.opt = True

        # collect and print the solution
        selected_cycles_chains = []
        for i in range(len(allo.cyclechains)):
            if allo.isActivated[i] == 1 and isCycleUsed[i].X > EPSILON:
                selected_cycles_chains.append(allo.cyclechains[i])

        return obj_val, selected_cycles_chains

    except gp.GurobiError as e:
        print(f"Gurobi Error: {e}")
        return -1, []
    except Exception as e:
        print(f"Exception during optimization: {e}")
        return -1, []
