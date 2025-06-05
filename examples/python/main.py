# examples/python_example/main.py
import sys
import os

# Add the build directory of the Python module to sys.path
# This is so Python can find the compiled .so/.pyd file.
# Adjust path based on your CMake build output structure.
# If _eos_fortran_legacy_lib.so is in build/python/:
module_build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../build/python'))
# If _eos_fortran_legacy_lib.so is directly in build/
# module_build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../build'))
print(f"Looking for Python module in: {module_build_dir}")
sys.path.insert(0, module_build_dir)

try:
    import PyEOS as eos_lib
    print("Successfully imported _eos_fortran_legacy_lib")
except ImportError as e:
    print(f"Error importing _eos_fortran_legacy_lib: {e}")
    print("Please ensure the module is built and its directory is in sys.path.")
    sys.exit(1)

def main():
    print("Python Example: Starting EOS Test")
    print("==================================")

    eos = eos_lib.EquationOfState()
    print("EquationOfState object created.")

    # Adjust data_dir path relative to this script or an absolute path
    # If script is in project_root/examples/python_example
    # and data is in project_root/eos_data_dir
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../eos_data_dir'))
    print(f"Using eos_data_dir: {data_dir}")

    eos_ids_to_load = [1, 2, 10000]

    try:
        print(f"\nInitializing EOS with IDs: {eos_ids_to_load}")
        init_stat = eos.initialize(eos_ids_to_load, data_dir)
        print(f"eos.initialize status: {init_stat}") # Should be EOS_SUCCESS (0)

        print("\nSetting control variables...")
        eos.set_complicated_eos_use_old_cold_term(True)
        eos.set_use_tfd_data_ver1(True)

        print("\nPerforming computations:")
        rho, T = 1.2, 300.0
        print(f"  Computing for EOS ID 1 (Air) with rho={rho}, T={T}")
        results_air = eos.compute(eos_id=1, rho=rho, T=T)
        print(f"    Results Air: P={results_air['P']:.5f}, E={results_air['E']:.5f}")

        rho, T = 1.8, 1200.0
        print(f"  Computing for EOS ID 10000 (Complicated) with rho={rho}, T={T}")
        results_comp = eos.compute(eos_id=10000, rho=rho, T=T)
        print(f"    Results Complicated: P={results_comp['P']:.5f}, E={results_comp['E']:.5f}")

        # Test pack and unpack
        print("\nTesting pack_data...")
        packed_bytes = eos.pack_data()
        print(f"  Packed data size: {len(packed_bytes)} bytes")
        if not packed_bytes:
            print("  Warning: pack_data returned empty bytes.")

        if packed_bytes:
            print("\nTesting unpack_data...")
            eos_new = eos_lib.EquationOfState() # New instance
            unpack_stat = eos_new.unpack_data(packed_bytes)
            # unpack_data in bindings currently returns istat or raises exception
            # print(f"  unpack_data status: {unpack_stat}")

            print(f"  Re-computing for ID 10000 on unpacked instance with rho={rho}, T={T}")
            results_new_comp = eos_new.compute(eos_id=10000, rho=rho, T=T)
            print(f"    Results New Complicated: P={results_new_comp['P']:.5f}, E={results_new_comp['E']:.5f}")

            # Compare P from original and new
            if abs(results_comp['P'] - results_new_comp['P']) < 1e-9:
                print("  Pack/Unpack Verified: P values match.")
            else:
                print("  ERROR: Pack/Unpack Verification FAILED: P values differ.")


    except RuntimeError as e:
        print(f"\nRuntimeError during EOS operations: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
    finally:
        if 'eos' in locals() and eos is not None:
            print("\nFreeing resources...")
            eos.free_resources() # Explicitly call if needed, or rely on destructor

    print("\n==================================")
    print("Python Example: Test Complete.")

if __name__ == "__main__":
    main()
