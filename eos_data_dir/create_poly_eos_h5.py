# create_poly_eos_h5.py
import h5py
import numpy as np

filename = "mat102/eos10200.h5"

with h5py.File(filename, 'w') as f:
    # Model Info
    model_info_grp = f.create_group("model_info")
    dt = h5py.string_dtype(encoding='utf-8')
    model_info_grp.create_dataset("type_name", data="polynomial_v1_hdf5", dtype=dt)
    model_info_grp.create_dataset("eos_id", data=np.int32(10200))
    model_info_grp.create_dataset("eos_name", data="Sample Poly EOS for 10200", dtype=dt)

    # Parameters
    poly_grp = f.create_group("polynomial")
    coeffs = np.array([1.0, 2.5, -0.3, 0.01], dtype=np.float64) # Example coefficients
    poly_grp.create_dataset("coefficients", data=coeffs)

    validity_grp = f.create_group("validity")
    validity_grp.create_dataset("rho_min", data=0.1)
    validity_grp.create_dataset("rho_max", data=10.0)
    validity_grp.create_dataset("temp_min", data=100.0)
    validity_grp.create_dataset("temp_max", data=5000.0)
print(f"Created {filename}")
