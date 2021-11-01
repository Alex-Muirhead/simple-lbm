#include "hdf5.h"
#include "lbm.h"
#include "hdfio.h"

// Function modified from Tutorial 1
int mat2hdf5(double *wdata, int shape_x, int shape_y)
{
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

class File {
   public:
    hid_t file;

    File open(const char *filename) {
        File newfile = File();
        newfile.file = H5Fopen(
            filename,
            H5F_ACC_RDONLY,  // Read Only access
            H5P_DEFAULT      // Default access properties
        );

        return newfile;
    }

    File create(const char *filename) {
        File newfile = File();
        newfile.file = H5Fcreate(
            filename,
            H5F_ACC_TRUNC,  // Truncate existing file if any
            H5P_DEFAULT,    // Default create properties
            H5P_DEFAULT     // Default access properties
        );

        return newfile;
    }

    void close() {
        // Close file resources
        H5Fclose(file);
    }

    ~File() {
        // If I forget
        close();
    }

    void read_double(const char *dataname, void *buffer) {
        hid_t dataset = H5Dopen(
            this->file,
            dataname,
            H5P_DEFAULT  // Unknown default options
        );

        H5Dclose(dataset);
    }
};