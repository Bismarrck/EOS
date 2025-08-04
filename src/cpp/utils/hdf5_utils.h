#ifndef HDF5_UTILS_H_H
#define HDF5_UTILS_H_H

#include <hdf5.h>

#include <string>
#include <vector>

namespace HDF5Utils {

// --- HDF5 Helper Implementations ---
herr_t read_hdf5_scalar_int(hid_t group_id_or_file_id, const char *dset_name,
                            int &out_val);

herr_t read_hdf5_scalar_double(hid_t group_id_or_file_id, const char *dset_name,
                               double &out_val);

herr_t read_hdf5_dataset_double(hid_t file_id, const char *dset_name,
                                int expected_n1, int expected_n2,
                                std::vector<double> &out_data);

// Helper for reading string attributes
bool read_hdf5_string_attribute(hid_t object_id, const char *attr_name,
                                std::string &out_str);

// Helper for reading scalar variable-length string dataset
bool read_hdf5_scalar_vlen_string_dataset(hid_t file_id, const char *dset_name,
                                          std::string &out_str);

herr_t write_hdf5_dataset_double_1d(hid_t file_id, const char *dset_name,
                                    const hsize_t *dims, const double *data);
herr_t write_hdf5_dataset_double_2d(hid_t file_id, const char *dset_name,
                                    const hsize_t *dims, const double *data);

}  // namespace HDF5Utils

#endif  // HDF5_UTILS_H_H
