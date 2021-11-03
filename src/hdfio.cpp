#include "hdfio.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "constants.h"
#include "hdf5.h"
#include "index.h"

#define CHECK_STATUS(x)                \
    if (x) {                           \
        std::cerr << "HDF5 Error in [" \
                  << __FILE__ << ":"   \
                  << __LINE__ << "]"   \
                  << endl;             \
    }

using namespace std;

// Function modified from Tutorial 1
int mat2hdf5(double *wdata, int shape_x, int shape_y) {
    hid_t file, space, dset; /* Handles */
    herr_t status;
    hsize_t dims[3] = {(hsize_t)shape_y, (hsize_t)shape_x, (hsize_t)Q};

    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate(WRITE_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple(3, dims, NULL);

    /*
     * Create the dataset and write the floating point data to it.  In
     * this example we will save the data as 64 bit little endian IEEE
     * floating point numbers, regardless of the native type.  The HDF5
     * library automatically converts between different floating point
     * types.
     */
    dset = H5Dcreate(file, DATASET, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      wdata);

    /*
     * Close and release resources.
     */
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);

    return 0;
}

/* ------------------- Input Data Structure ------------------- */

InputData::InputData(hid_t file) {
    hid_t status, attribute;

    this->file = file;

    printf("Loading configuration data\n");

    attribute = H5Aopen_by_name(
        file,
        "layout",
        "size_x",
        H5P_DEFAULT,
        H5P_DEFAULT);

    status = H5Aread(attribute, H5T_NATIVE_INT, &shape_x);
    CHECK_STATUS(status)
    H5Aclose(attribute);

    printf("Loaded size_x=%d\n", shape_x);

    attribute = H5Aopen_by_name(
        file,
        "layout",
        "size_y",
        H5P_DEFAULT,
        H5P_DEFAULT);

    status = H5Aread(attribute, H5T_NATIVE_INT, &shape_y);
    CHECK_STATUS(status)
    H5Aclose(attribute);

    printf("Loaded size_y=%d\n", shape_y);
}

InputData InputData::open(const char *filename) {
    hid_t file = H5Fopen(
        filename,
        H5F_ACC_RDONLY,  // Read Only access
        H5P_DEFAULT      // Default access properties
    );

    return InputData(file);
}

void InputData::close() {
    // Close file resources
    H5Fclose(file);
}

void InputData::read_scalar(const char *dataname, double *buffer) {
    hid_t status, dataset;

    dataset = H5Dopen(
        file,
        dataname,
        H5P_DEFAULT  // Default access properties
    );

    status = H5Dread(
        dataset,
        H5T_NATIVE_DOUBLE,  // Datatype
        H5S_ALL,            // Write to all of memory buffer
        H5S_ALL,            // Read in all of dataset / implicit dataspace
        H5P_DEFAULT,        // Default read properties
        buffer);
    CHECK_STATUS(status)

    // Clear up resources
    H5Dclose(dataset);
}

void InputData::configure_domain(int* buffer) {
    hid_t status, dataset;

    dataset = H5Dopen(
        file,
        "layout/types",
        H5P_DEFAULT  // Default access properties
    );

    status = H5Dread(
        dataset,
        H5T_NATIVE_INT,  // Datatype
        H5S_ALL,            // Write to all of memory buffer
        H5S_ALL,            // Read in all of dataset / implicit dataspace
        H5P_DEFAULT,        // Default read properties
        buffer);
    CHECK_STATUS(status)

    // Clear up resources
    H5Dclose(dataset);
}

/* ------------------- Output Data Structure ------------------- */

OutputData OutputData::open(const char *filename) {
    OutputData newfile = OutputData();
    newfile.file = H5Fopen(
        filename,
        H5F_ACC_RDONLY,  // Read Only access
        H5P_DEFAULT      // Default access properties
    );

    return newfile;
}

OutputData OutputData::create(const char *filename) {
    OutputData newfile = OutputData();
    newfile.file = H5Fcreate(
        filename,
        H5F_ACC_TRUNC,  // Truncate existing file if any
        H5P_DEFAULT,    // Default create properties
        H5P_DEFAULT     // Default access properties
    );

    return newfile;
}

void OutputData::close() {
    // Close file resources
    for (auto const &pair : datasets) {
        H5Dclose(pair.second);
    }

    H5Fclose(file);
}

// For initialising datasets
void OutputData::write_scalar(const char *dataname, double *data) {
    hid_t status, dataset, dataspace;

    // Temporary values
    unsigned int size_x = 128;
    unsigned int size_y = 128;

    int rank = 3;
    hsize_t dims[3] = {0, size_x, size_y};
    hsize_t maxdims[3] = {H5S_UNLIMITED, size_x, size_y};

    dataspace = H5Screate_simple(
        rank,
        dims,
        maxdims);

    hsize_t chunk_dims[3] = {1, size_x, size_y};  // Will only change the number of slices
    hid_t chunk_params = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(
        chunk_params,
        rank,
        chunk_dims);
    CHECK_STATUS(status)

    dataset = H5Dcreate(
        file,
        dataname,
        H5T_NATIVE_DOUBLE,
        dataspace,
        H5P_DEFAULT,
        chunk_params,
        H5P_DEFAULT);

    H5Pclose(chunk_params);

    // Save the dataset id to the map
    datasets.insert(pair<string, hid_t>(dataname, dataset));

    // Use existing function to add data
    append_scalar(dataname, data);
}

void OutputData::append_scalar(const char *dataname, double *data) {
    hid_t status, dataset, dataspace, memspace;
    hsize_t memdims[2], dims[3], maxdims[3], offset[3], count[3];

    // Temporary values
    unsigned int size_x = 128;
    unsigned int size_y = 128;

    // Get the relevant dataset
    dataset = datasets.at(dataname);
    dataspace = H5Dget_space(dataset);

    H5Sget_simple_extent_dims(dataspace, dims, maxdims);

    dims[0] += 1;  // Increase dataset size by 1 slice/chunk
    status = H5Dextend(dataset, dims);
    CHECK_STATUS(status)

    memdims[0] = size_x;
    memdims[1] = size_y;
    memspace = H5Screate_simple(2, memdims, NULL);

    dataspace = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace, dims, maxdims);

    offset[0] = dims[0] - 1;  // Whatever we're currently up to
    offset[1] = 0;
    offset[2] = 0;

    count[0] = 1;  // Only select a single slice
    count[1] = size_x;
    count[2] = size_y;

    status = H5Sselect_hyperslab(
        dataspace,
        H5S_SELECT_SET,  // Set the dataspace selection to this
        offset,          // The offset from origin of dataset
        NULL,            // Stride = 1 (Default)
        count,           // The number of data points in each direction
        NULL             // Blocksize = 1, no blocks (Default)
    );
    CHECK_STATUS(status)

    status = H5Dwrite(
        dataset,
        H5T_NATIVE_DOUBLE,  // Datatype to write
        memspace,           // Space within the memory buffer
        dataspace,          // Space within the dataset
        H5P_DEFAULT,        // Default property list
        data);
    CHECK_STATUS(status)

    H5Sclose(memspace);
    H5Sclose(dataspace);
}