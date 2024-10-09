import sys
import os
from allocation import Allocation as NormalAllocation
from allocation_generalized import Allocation as GeneralizedAllocation
from cycle_chain_deactivation import run_cycle_chain_deactivation as normal_run
from cycle_chain_deactivation_generalized import (
    run_cycle_chain_deactivation as generalized_run,
)


def import_kidney_data(filepath):
    """
    Imports kidney exchange data from a text file and structures it into a dictionary.

    Args:
        filepath (str): The path to the text file containing the kidney exchange data.

    Returns:
        dict: A dictionary containing the following keys:
            - num_pairs (int): The number of pairs (donor-patient pairs) in the instance.
            - num_ndd (int): The number of non-directed donors (NDDs) in the instance.
            - num_arcs (int): The total number of arcs in the instance.
            - pairs (list): A list of dictionaries, each representing a pair or NDD with the following keys:
                - id (int): The unique ID of the pair.
                - is_ndd (bool): True if the pair is an NDD, False otherwise.
                - donor_blood_type (int): Donor's blood type (0 = A, 1 = B, 2 = AB, 3 = O).
                - patient_blood_type (int): Patient's blood type (0 = A, 1 = B, 2 = AB, 3 = O).
                - patient_vpra (int): Patient's vPRA score (0 = below 0.5, 1 = between 0.5 and 0.85, 2 = above 0.85).
            - arcs (list): A list of dictionaries representing the arcs between pairs with the following keys:
                - donor_id (int): ID of the donor pair.
                - patient_id (int): ID of the patient pair.
                - weight (int): The weight of the arc (always 1 in this case).

    Example:
        data = import_kidney_data('path/to/data.txt')
        print(data['num_pairs'])  # Output: Number of pairs in the instance
    """

    data = {}

    # Open the file and read all lines into memory
    with open(filepath, "r") as f:
        lines = f.readlines()

        # Parse metadata: number of pairs, NDDs, and arcs
        num_pairs = int(lines[0].split(" ")[2])
        num_ndd = int(lines[1].split(" ")[2])
        num_arcs = int(lines[2].split(" ")[2])

        # Store the parsed values in the data dictionary
        data["num_pairs"] = num_pairs
        data["num_ndd"] = num_ndd
        data["num_arcs"] = num_arcs

        # Total number of entities (pairs + NDDs)
        num_things = num_pairs + num_ndd

        # List to store information about pairs and NDDs
        pairs = []

        # Parse each line corresponding to pairs and NDDs
        for line in lines[3 : num_things + 3]:
            id, is_ndd, donor_blood_type, patient_blood_type, patient_vpra = map(
                int, line.strip().split(",")
            )

            # Add pair/NDD info to the list
            pairs.append(
                {
                    "id": id,
                    "is_ndd": bool(is_ndd),
                    "donor_blood_type": donor_blood_type,
                    "patient_blood_type": patient_blood_type,
                    "patient_vpra": patient_vpra,
                }
            )

        # Store the pairs data in the dictionary
        data["pairs"] = pairs

        # List to store information about arcs
        arcs = []

        # Parse each line corresponding to arcs (donor-patient relationships)
        for line in lines[num_things + 3 :]:
            arc, weight = line.strip().split(",1,")

            # Extract donor and patient IDs from the arc
            donor_id, patient_id = arc.split(",")
            weight = int(weight.strip())  # Weight is always 1 in this context
            donor_id = int(donor_id[1:])  # Remove prefix and convert to int
            patient_id = int(patient_id[:-1])  # Remove suffix and convert to int

            # Add arc info to the list
            arcs.append(
                {
                    "donor_id": donor_id,
                    "patient_id": patient_id,
                    "weight": weight,
                }
            )

        # Store the arcs data in the dictionary
        data["arcs"] = arcs

    # Return the final structured dictionary
    return data


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Usage: python3 run.py <input_file> <output_file> [-g] [-c <max_cycle_length>] [-h <max_chain_length>]"
        )
        sys.exit(1)

    file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    # Set default values for cycle and chain lengths
    generalized = False
    max_cycle_length = 3
    max_chain_length = 4

    # Parse additional arguments for generalized version
    if "-g" in sys.argv:
        generalized = True
        if "-c" in sys.argv:
            max_cycle_length = int(sys.argv[sys.argv.index("-c") + 1])
        if "-h" in sys.argv:
            max_chain_length = int(sys.argv[sys.argv.index("-h") + 1])

    data = import_kidney_data(file_path)

    if generalized:
        print(
            f"Running generalized version with max cycle length {max_cycle_length} and max chain length {max_chain_length}..."
        )
        generalized_run(
            data,
            output_file_path,
            max_cycle_length=max_cycle_length,
            max_chain_length=max_chain_length,
        )
    else:
        print(
            "Running normal version with fixed max cycle length 3 and max chain length 4..."
        )
        normal_run(data, output_file_path)
