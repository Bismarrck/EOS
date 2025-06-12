# In create_sample_tfd.py
import h5py
import os
import numpy as np
import datetime

def create_tfd_h5_with_metadata(filename, n1, n2, matrix_a_data, matrix_b_data,
                                model_desc, matrix_a_desc, matrix_b_desc):
    if os.path.exists(filename):
        os.remove(filename)

    with h5py.File(filename, 'w') as f:
        # File level attributes
        f.attrs['file_creation_date'] = str(datetime.date.today())
        f.attrs['tfd_model_version'] = "1.1" # Example version

        # Common dimensions
        common_dims_grp = f.create_group("common_dimensions")
        common_dims_grp.create_dataset("N1", data=np.int32(n1))
        common_dims_grp.create_dataset("N2", data=np.int32(n2))

        # Model Description (Scalar Variable-Length String Dataset)
        # For variable length unicode strings, h5py needs special dtype
        dt = h5py.string_dtype(encoding='utf-8') # Use 'ascii' if only ASCII
        f.create_dataset("model_description", data=model_desc, dtype=dt)

        # Matrix A
        dset_a = f.create_dataset("matrix_A", data=matrix_a_data, dtype='f8', compression="gzip")
        dset_a.attrs['description'] = matrix_a_desc
        dset_a.attrs['units'] = "kJ/mol" # Example

        # Matrix B
        dset_b = f.create_dataset("matrix_B", data=matrix_b_data, dtype='f8', compression="gzip")
        dset_b.attrs['description'] = matrix_b_desc

        print(f"Created HDF5 TFD file with metadata: {filename}")

if __name__ == "__main__":
    N1_val, N2_val = 100, 50
    np.random.seed(1)

    model_description_text = (
            "This TFD file contains matrices for a sample material.\n"
            "Version 1.1 includes updated coefficients for high temperature.\n"
            "Generated on: " + str(datetime.date.today())
    )
    mat_a_desc_text = "Primary energy contribution matrix A."
    mat_b_desc_text = "Entropy contribution matrix A."

    # Sample data for tfd_ver1.h5
    mat_a_v1 = np.random.rand(N1_val, N2_val)
    mat_a_v1[0, 0: 6] = 1.1, 1.2, 1.3, 2.1, 2.2, 2.3
    mat_b_v1 = np.random.rand(N1_val, N2_val)
    mat_b_v1[0, 0: 6] = 10.1, 10.2, 10.3, 20.1, 20.2, 20.3
    create_tfd_h5_with_metadata("tfd_ver1.h5", N1_val, N2_val, mat_a_v1, mat_b_v1,
                                model_description_text, mat_a_desc_text, mat_b_desc_text)

    model_description_text = (
            "This TFD file contains matrices for a sample material.\n"
            "Version 1.1 includes updated coefficients for high temperature.\n"
            "Generated on: " + str(datetime.date.today())
    )
    mat_a_desc_text = "Primary energy contribution matrix B."
    mat_b_desc_text = "Entropy contribution matrix B."

    mat_a_v2 = np.random.rand(N1_val, N2_val) * 15.0 + 5.0
    mat_b_v2 = np.random.rand(N1_val, N2_val) * 25.0 + 3.0
    create_tfd_h5_with_metadata("tfd_ver2.h5", N1_val, N2_val, mat_a_v2, mat_b_v2,
                                model_description_text, mat_a_desc_text, mat_b_desc_text)
