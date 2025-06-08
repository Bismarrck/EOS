# examples/python_example/test_pyeos.py
import sys
import os
import platform

# --- Helper to find the compiled PyEOS module ---

# Determine the build directory path relative to this script's location
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

print(f"Successfully imported PyEOS version (pybind11 version): {PyEOS.__doc__}") # __doc__ from module

def print_py_eos_results(eos_id, results_tuple, istat_implicit=0):
    # results_tuple is (P, E, dPdT, dEdT, dPdrho)
    # istat is handled by exceptions in the binding layer
    print(f"  EOS ID {eos_id} SUCCESS:")
    print(f"    P      = {results_tuple[0]:.5f}")
    print(f"    E      = {results_tuple[1]:.5f}")
    print(f"    dPdT   = {results_tuple[2]:.5f}")
    print(f"    dEdT   = {results_tuple[3]:.5f}")
    print(f"    dPdrho = {results_tuple[4]:.5f}")

def main():
    print("\nPython Example: Starting EOS Test")
    print("==================================")

    eos_manager = PyEOS.EquationOfState()

    # Define eos_data_dir. Adjust path as necessary.
    # Path relative to this script, assuming script is in examples/python_example
    # and eos_data_dir is at project root.
    eos_data_path = os.path.join(project_root, "eos_data_dir")
    print(f"Using data_dir: {eos_data_path}")

    # --- Check Data Directory ---
    try:
        eos_manager.check_eos_data_dir(eos_data_path, [10000]) # Check one complicated ID
        print("Python: eos_check_data_dir successful.")
    except RuntimeError as e:
        print(f"Python Error during check_eos_data_dir: {e}")
        return

    # --- Initialize EOS Instance ---
    eos_ids_to_load = [1, 2, 10000]
    try:
        eos_manager.initialize(eos_ids_to_load, eos_data_path)
        print(f"Python: EOS Initialized with IDs {eos_ids_to_load}.")
    except RuntimeError as e:
        print(f"Python Error during initialize: {e}")
        return

    # --- Set a Control Variable ---
    try:
        eos_manager.set_complicated_eos_use_old_cold_term(True)
        print("Python: set_complicated_eos_use_old_cold_term(True) successful.")
        eos_manager.set_use_tfd_data_ver1(True) # Explicitly TFD v1
        print("Python: set_use_tfd_data_ver1(True) successful.")
    except RuntimeError as e:
        print(f"Python Error setting control variable: {e}")
        return

    # --- Perform Computations ---
    print("\nPerforming Computations:")
    try:
        # Analytic Air (ID 1)
        rho, temp = 1.2, 300.0
        print(f"Computing for EOS ID 1 (Air Analytic) with rho={rho}, T={temp}")
        results = eos_manager.compute(1, rho, temp)
        print_py_eos_results(1, results)

        # Complicated EOS (ID 10000)
        rho, temp = 1.8, 1200.0
        print(f"Computing for EOS ID 10000 (Complicated) with rho={rho}, T={temp}")
        results = eos_manager.compute(10000, rho, temp)
        print_py_eos_results(10000, results)

        # Test error case for Analytic Carbon (ID 2)
        rho, temp = -1.0, 500.0
        print(f"Computing for EOS ID 2 (Carbon Analytic) with rho={rho}, T={temp} (expecting error)")
        try:
            results = eos_manager.compute(2, rho, temp)
            print_py_eos_results(2, results) # Should not reach here
        except RuntimeError as e_compute:
            print(f"  Python: Caught expected error for ID 2: {e_compute}")

    except RuntimeError as e:
        print(f"Python Error during compute: {e}")
        return

    # --- Pack and Unpack Data ---
    print("\nPacking and Unpacking Data:")
    try:
        print("Packing data...")
        packed_bytes = eos_manager.pack_data()
        print(f"Packed data size: {len(packed_bytes)} bytes")

        print("Unpacking data into a new manager instance...")
        eos_manager_rankN = PyEOS.EquationOfState()
        eos_manager_rankN.unpack_data(packed_bytes)
        print("Unpack successful. Verifying by re-computing ID 10000:")

        rho, temp = 1.8, 1200.0 # Same inputs as before pack
        results_original = eos_manager.compute(10000, rho, temp) # Compute on original
        results_unpacked = eos_manager_rankN.compute(10000, rho, temp) # Compute on unpacked

        print("Original manager computation for ID 10000:")
        print_py_eos_results(10000, results_original)
        print("Unpacked manager computation for ID 10000:")
        print_py_eos_results(10000, results_unpacked)

        # Basic check if P and E match (use a tolerance for floats)
        tolerance = 1e-9
        if (abs(results_original[0] - results_unpacked[0]) < tolerance and
                abs(results_original[1] - results_unpacked[1]) < tolerance):
            print("Verification: Pack/Unpack successful, computation results match for ID 10000.")
        else:
            print("Verification ERROR: Pack/Unpack results DO NOT match for ID 10000.")


    except RuntimeError as e:
        print(f"Python Error during pack/unpack: {e}")
        return

    # --- Free Resources (optional in Python due to RAII in C++, but good if explicit free is there) ---
    try:
        eos_manager.free_resources()
        print("\nPython: Original EOS manager resources freed.")
        eos_manager_rankN.free_resources()
        print("Python: Unpacked EOS manager resources freed.")
    except RuntimeError as e:
        print(f"Python Error during free_resources: {e}")


    print("\n==================================")
    print("Python Example: Test Complete.")

if __name__ == "__main__":
    main()