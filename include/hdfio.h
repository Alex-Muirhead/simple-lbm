#define READ_FILE "lattice_in.h5"
#define WRITE_FILE "lattice_out.h5"
#define DATASET "DS1"

int mat2hdf5(double *wdata, int shape_x, int shape_y);