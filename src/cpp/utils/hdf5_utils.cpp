#include "hdf5_utils.h"

#include <hdf5.h>

#include <iostream>

// --- HDF5 Helper Implementations ---
herr_t HDF5Utils::read_hdf5_scalar_int(hid_t group_id_or_file_id,
                                       const char *dset_name, int &out_val) {
  hid_t dset_id = H5Dopen2(group_id_or_file_id, dset_name, H5P_DEFAULT);
  if (dset_id < 0) {
    std::cerr << "HDF5 Error: Could not open integer dataset '" << dset_name
              << "'" << std::endl;
    return -1;  // Indicate error
  }
  herr_t status =
      H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &out_val);
  if (status < 0) {
    std::cerr << "HDF5 Error: Could not read integer dataset '" << dset_name
              << "'" << std::endl;
  }
  H5Dclose(dset_id);
  return status;
}

herr_t HDF5Utils::read_hdf5_dataset_double(hid_t file_id, const char *dset_name,
                                           int expected_n1, int expected_n2,
                                           std::vector<double> &out_data) {
  hid_t dset_id = H5Dopen2(file_id, dset_name, H5P_DEFAULT);
  if (dset_id < 0) {
    std::cerr << "HDF5 Error: Could not open double dataset '" << dset_name
              << "'" << std::endl;
    return -1;
  }

  hid_t dataspace_id = H5Dget_space(dset_id);
  if (dataspace_id < 0) {
    std::cerr << "HDF5 Error: Could not get dataspace for '" << dset_name << "'"
              << std::endl;
    H5Dclose(dset_id);
    return -1;
  }

  int rank = H5Sget_simple_extent_ndims(dataspace_id);
  if (rank != 2) {  // Assuming 2D matrices
    std::cerr << "HDF5 Error: Dataset '" << dset_name
              << "' is not 2D (rank=" << rank << ")" << std::endl;
    H5Sclose(dataspace_id);
    H5Dclose(dset_id);
    return -1;
  }

  hsize_t dims[2];
  H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
  if (static_cast<int>(dims[0]) != expected_n1 ||
      static_cast<int>(dims[1]) != expected_n2) {
    std::cerr << "HDF5 Error: Dimension mismatch for dataset '" << dset_name
              << "'. "
              << "Expected " << expected_n1 << "x" << expected_n2 << ", Got "
              << dims[0] << "x" << dims[1] << std::endl;
    H5Sclose(dataspace_id);
    H5Dclose(dset_id);
    return -1;
  }

  out_data.resize(static_cast<size_t>(expected_n1) * expected_n2);
  herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, out_data.data());
  if (status < 0) {
    std::cerr << "HDF5 Error: Could not read double dataset '" << dset_name
              << "'" << std::endl;
  }

  H5Sclose(dataspace_id);
  H5Dclose(dset_id);
  return status;
}

// Helper for reading string attributes
bool HDF5Utils::read_hdf5_string_attribute(hid_t object_id,
                                           const char *attr_name,
                                           std::string &out_str) {
  out_str.clear();
  if (H5Aexists(object_id, attr_name) <= 0) {
    return true;  // Attribute doesn't exist, not an error for this helper
  }

  hid_t attr_id = H5Aopen_name(object_id, attr_name);
  if (attr_id < 0) { /* error */
    return false;
  }

  hid_t file_type_id = H5Aget_type(attr_id);
  if (file_type_id < 0) { /* error */
    H5Aclose(attr_id);
    return false;
  }

  // Check if it's a variable-length string in the file
  htri_t is_vlen_str_in_file = H5Tis_variable_str(file_type_id);
  if (is_vlen_str_in_file < 0) { /* error */
    H5Tclose(file_type_id);
    H5Aclose(attr_id);
    return false;
  }

  hid_t mem_type_id = H5Tcopy(H5T_C_S1);  // Start with a C string type
  if (mem_type_id < 0) {                  /* error */
    H5Tclose(file_type_id);
    H5Aclose(attr_id);
    return false;
  }

  if (is_vlen_str_in_file > 0) {  // Variable-length string in file
    if (H5Tset_size(mem_type_id, H5T_VARIABLE) < 0) { /* error */
      H5Tclose(mem_type_id);
      H5Tclose(file_type_id);
      H5Aclose(attr_id);
      return false;
    }
  } else if (H5Tget_class(file_type_id) ==
             H5T_STRING) {  // Fixed-length string in file
    size_t fixed_size = H5Tget_size(file_type_id);
    if (fixed_size == 0 || fixed_size == H5T_VARIABLE) { /* error: unexpected */
      H5Tclose(mem_type_id);
      H5Tclose(file_type_id);
      H5Aclose(attr_id);
      return false;
    }
    if (H5Tset_size(mem_type_id, fixed_size) < 0) { /* error */
      H5Tclose(mem_type_id);
      H5Tclose(file_type_id);
      H5Aclose(attr_id);
      return false;
    }
  } else { /* error: not a string type */
    H5Tclose(mem_type_id);
    H5Tclose(file_type_id);
    H5Aclose(attr_id);
    return false;
  }

  // *** Crucial: Match character set if reading UTF-8 from file ***
  // h5dump shows CSET H5T_CSET_UTF8 for the attributes.
  if (H5Tget_cset(file_type_id) == H5T_CSET_UTF8) {
    if (H5Tset_cset(mem_type_id, H5T_CSET_UTF8) < 0) {
      std::cerr << "HDF5 Warning: Failed to set memory type character set to "
                   "UTF-8 for attribute '"
                << attr_name << "'." << std::endl;
      // Continue anyway, but conversion might still fail or be lossy.
    }
  }
  // Ensure null termination for memory type (good practice)
  // H5Tset_strpad(mem_type_id, H5T_STR_NULLTERM); // Usually H5T_C_S1 is
  // already null-terminated

  char *str_read_buffer_vlen = nullptr;
  std::vector<char> str_read_buffer_fixed;
  void *actual_read_ptr;

  if (is_vlen_str_in_file > 0) {
    actual_read_ptr = &str_read_buffer_vlen;  // H5Aread expects char**
  } else {
    size_t mem_size = H5Tget_size(mem_type_id);
    str_read_buffer_fixed.resize(mem_size + 1,
                                 '\0');  // +1 for safety null term
    actual_read_ptr = str_read_buffer_fixed.data();
  }

  herr_t read_status = H5Aread(attr_id, mem_type_id, actual_read_ptr);

  if (read_status < 0) {
    std::cerr << "HDF5 Error: H5Aread failed for attribute '" << attr_name
              << "'." << std::endl;
    H5Eprint(H5E_DEFAULT, stderr);
    if (is_vlen_str_in_file > 0 && str_read_buffer_vlen)
      H5free_memory(str_read_buffer_vlen);
    H5Tclose(mem_type_id);
    H5Tclose(file_type_id);
    H5Aclose(attr_id);
    return false;
  }

  if (is_vlen_str_in_file > 0) {
    if (str_read_buffer_vlen) {
      out_str = str_read_buffer_vlen;
      H5free_memory(str_read_buffer_vlen);
    } else {
      out_str = "";
    }
  } else {
    out_str = str_read_buffer_fixed.data();
  }

  H5Tclose(mem_type_id);
  H5Tclose(file_type_id);
  H5Aclose(attr_id);
  return true;
}

// Helper for reading scalar variable-length string dataset
bool HDF5Utils::read_hdf5_scalar_vlen_string_dataset(hid_t file_id,
                                                     const char *dset_name,
                                                     std::string &out_str) {
  out_str.clear();
  hid_t dset_id = H5Dopen2(file_id, dset_name, H5P_DEFAULT);
  if (dset_id < 0) { /* ... */
    return false;
  }

  hid_t file_type_id = H5Dget_type(dset_id);
  if (file_type_id < 0) { /* ... */
    H5Dclose(dset_id);
    return false;
  }

  // Verify it's a variable-length string in the file
  htri_t is_vlen_str_in_file = H5Tis_variable_str(file_type_id);
  if (is_vlen_str_in_file <= 0) { /* ... error message ... */
    H5Tclose(file_type_id);
    H5Dclose(dset_id);
    return false;
  }

  // Verify it's scalar
  hid_t dspace_id = H5Dget_space(dset_id);
  if (dspace_id < 0 ||
      H5Sget_simple_extent_type(dspace_id) != H5S_SCALAR) { /* ... */
    if (dspace_id >= 0) H5Sclose(dspace_id);
    H5Tclose(file_type_id);
    H5Dclose(dset_id);
    return false;
  }
  H5Sclose(dspace_id);

  // Prepare memory type
  hid_t mem_type_id = H5Tcopy(H5T_C_S1);
  if (mem_type_id < 0) { /* error */
    H5Tclose(file_type_id);
    H5Dclose(dset_id);
    return false;
  }
  if (H5Tset_size(mem_type_id, H5T_VARIABLE) < 0) { /* error */
    H5Tclose(mem_type_id);
    H5Tclose(file_type_id);
    H5Dclose(dset_id);
    return false;
  }

  // *** Match character set for dataset as well ***
  if (H5Tget_cset(file_type_id) ==
      H5T_CSET_UTF8) {  // From h5dump, dataset is also UTF8
    if (H5Tset_cset(mem_type_id, H5T_CSET_UTF8) < 0) {
      std::cerr << "HDF5 Warning: Failed to set memory type character set to "
                   "UTF-8 for dataset '"
                << dset_name << "'." << std::endl;
    }
  }

  char *str_read_buffer = nullptr;
  herr_t read_status = H5Dread(dset_id, mem_type_id, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, &str_read_buffer);

  if (read_status < 0) {
    std::cerr << "HDF5 Error: H5Dread failed for dataset '" << dset_name << "'."
              << std::endl;
    H5Eprint(H5E_DEFAULT, stderr);
    if (str_read_buffer) H5free_memory(str_read_buffer);
    H5Tclose(mem_type_id);
    H5Tclose(file_type_id);
    H5Dclose(dset_id);
    return false;
  }

  if (str_read_buffer) {
    out_str = str_read_buffer;
    H5free_memory(str_read_buffer);
  } else {
    out_str = "";
  }

  H5Tclose(mem_type_id);
  H5Tclose(file_type_id);
  H5Dclose(dset_id);
  return true;
}
