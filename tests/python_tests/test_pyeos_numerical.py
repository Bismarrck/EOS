import pytest
import csv
import math
import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, "..", "..")) # Up two levels to project root

# Common CMake build directory names
build_dir_names = ["build", "cmake-build-debug", "cmake-build-release", "_build"]
# Suffix for python modules
module_dir_name = "python"

# Try to find the build path that contains the python module
build_path_found = None
for build_name in build_dir_names:
    potential_path = os.path.join(project_root, build_name, module_dir_name)
    if os.path.isdir(potential_path):
        build_path_found = potential_path
        break

if build_path_found:
    print(f"Found Python module build path: {build_path_found}")
    sys.path.insert(0, build_path_found)
else:
    print("Error: Could not automatically find the Python module build directory.")
    print("Please ensure the PyEOS module is built and its directory is in PYTHONPATH,")
    print(f"or run this script from a directory where 'import PyEOS' works.")
    print(f"Searched common build paths like '{os.path.join(project_root, 'build', module_dir_name)}'")
    sys.exit(1)

try:
    import PyEOS
except ImportError as e:
    print(f"Error importing PyEOS: {e}")
    print("Make sure the module is built and its location is in sys.path or PYTHONPATH.")
    sys.exit(1)

script_dir_py = os.path.dirname(os.path.abspath(__file__))
project_root_py = os.path.abspath(os.path.join(script_dir_py, "..", ".."))
REF_DATA_FILE = os.path.join(project_root_py, "tests", "test_data", "reference_eos_data.csv")
EOS_DATA_DIR_PY = os.path.join(project_root_py, "eos_data_dir")

def load_ref_cases_py(filepath):
    cases = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(filter(lambda row: row.strip() and not row.startswith('#'), f))
        for row in reader:
            try:
                # Basic type conversion, add error handling
                case = {
                    'eos_id': int(row['eos_id']),
                    'rho': float(row['rho'].replace('D','E').replace('d','e')),
                    'T': float(row['T'].replace('D','E').replace('d','e')),
                    'ref_P': float(row['ref_P'].replace('D','E').replace('d','e')),
                    'ref_E': float(row['ref_E'].replace('D','E').replace('d','e')),
                    'ref_dPdT': float(row['ref_dPdT'].replace('D','E').replace('d','e')),
                    'ref_dEdT': float(row['ref_dEdT'].replace('D','E').replace('d','e')),
                    'ref_dPdrho': float(row['ref_dPdrho'].replace('D','E').replace('d','e')),
                    'description': row.get('description', f"ID {row['eos_id']}"),
                    'use_tfd_v1': row.get('use_tfd_v1', 'true').lower() in ['true', '1', '']
                }
                cases.append(case)
            except Exception as e:
                print(f"Error parsing row in Python: {row} -> {e}")
    return cases

test_cases_py = load_ref_cases_py(REF_DATA_FILE)


@pytest.mark.parametrize("tc", test_cases_py, ids=lambda tc: tc['description'])
def test_eos_numerical_py(tc):
    eos = PyEOS.EquationOfState()

    eos.set_use_tfd_data_ver1(tc['use_tfd_v1'])
    # Initialize with only the current eos_id needed for this test case to speed things up
    # or initialize with all unique IDs once if that's preferred and feasible.
    try:
        eos.initialize([tc['eos_id']], EOS_DATA_DIR_PY)
    except RuntimeError as e:
        pytest.fail(f"Initialization failed for {tc['description']}: {e}")

    try:
        P_calc, E_calc, dPdT_calc, dEdT_calc, dPdrho_calc = eos.compute(tc['eos_id'], tc['rho'], tc['T'])
    except RuntimeError as e:
        pytest.fail(f"Compute failed for {tc['description']}: {e}")

    rel_tol = 1e-7
    abs_tol = 1e-9
    assert math.isclose(P_calc, tc['ref_P'], rel_tol=rel_tol, abs_tol=abs_tol), f"P mismatch for {tc['description']}"
    assert math.isclose(E_calc, tc['ref_E'], rel_tol=rel_tol, abs_tol=abs_tol), f"E mismatch for {tc['description']}"
