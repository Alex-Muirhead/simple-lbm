#include <array>
#include <iostream>
#include <sstream>

#include "hdf5.h"
#include "hdfio.h"
#include "lbm.h"

using namespace std;

class Field {
   private:
    static int shape_x;
    static int shape_y;
    static int shape_q;

   public:
    static int size;
    static void set_field_shape(int shape_x, int shape_y, int shape_q) {
        Field::shape_x = shape_x;
        Field::shape_y = shape_y;
        Field::shape_q = shape_q;
        Field::size = shape_x * shape_y * shape_q;
    }

    static unsigned int index(int x, int y, int q) {
        return shape_q * (shape_x * y + x) + q;
    }
};

int Field::shape_x;
int Field::shape_y;
int Field::shape_q;
int Field::size;

class Scalar {
   private:
    static int shape_x;
    static int shape_y;

   public:
    static int size;
    static void set_scalar_shape(int shape_x, int shape_y) {
        Scalar::shape_x = shape_x;
        Scalar::shape_y = shape_y;
        Scalar::size = shape_x * shape_y;
    }

    static unsigned int index(int x, int y) {
        return shape_x * y + x;
    }

    static double* create() {
        return new double[shape_x * shape_y];
    }
};

int Scalar::shape_x;
int Scalar::shape_y;
int Scalar::size;

/*
 * Ensures a positive value from modulo
 * Function taken from https://stackoverflow.com/a/23214219
 */
inline int mod(int k, int n) {
    return ((k %= n) < 0) ? k + n : k;
}

enum PointType : int {
    Fluid  = 1 << 0,   // Standard streaming
    NoSlip = 1 << 1,   // Bounce-back boundary
    Slip_V = 1 << 2,   // Vertical bounce-forward boundary
    Slip_H = 1 << 3,   // Horiztonal bounce-forward boundary
};

int main(int argc, char* argv[]) {
    int val;
    if (argc >= 2) {
        std::istringstream iss(argv[1]);

        if (iss >> val) {
            // Conversion successful
        }
    }

    double nu_lb = 0.092;                // Lattice dynamic viscosity
    double omega = 1.0 / (3. * nu_lb + 0.5);  // Relaxation parameter

    // HACKING IN VALUES

    hid_t file, dset, space; /* Handles */
    herr_t status;

    /*
     * Open an existing file using the default properties.
     */
    file = H5Fopen("taylor_green.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

    /*
     * Create the dataset and read the floating point data from it. The HDF5
     * library automatically converts between different floating point types.
     */
    dset = H5Dopen(file, "types", H5P_DEFAULT);

    hsize_t dims[2];

    space = H5Dget_space(dset);
    H5Sget_simple_extent_dims(space, dims, NULL);

    int size_y = dims[0];
    int size_x = dims[1];

    Field::set_field_shape(size_x, size_y, Q);
    Scalar::set_scalar_shape(size_x, size_y);

    int* type_lattice = new int[Scalar::size];
    status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, type_lattice);

    status = H5Dclose(dset);
    status = H5Sclose(space);

    double* rho   = new double[Scalar::size];
    double* vel_x = new double[Scalar::size];
    double* vel_y = new double[Scalar::size];

    dset = H5Dopen(file, "initial/density", H5P_DEFAULT);
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho);

    status = H5Dclose(dset);

    dset = H5Dopen(file, "initial/vel_x", H5P_DEFAULT);
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_x);

    status = H5Dclose(dset);

    dset = H5Dopen(file, "initial/vel_y", H5P_DEFAULT);
    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_y);

    status = H5Dclose(dset);
    status = H5Fclose(file);

    double* f1 = new double[Field::size];
    double* f2 = new double[Field::size];

    // Set up equilibrium values for initial flow
    for (int y = 0; y < size_y; y++) {
        for (int x = 0; x < size_x; x++) {
            // TODO: Use an imported file instead
            if (type_lattice[Scalar::index(x, y)] != PointType::Fluid) {
                // Set boundary lattice at zero
                for (int q = 0; q < Q; q++) {
                    f1[Field::index(x, y, q)] = 1.0;  // Initial density
                }
                continue;
            }

            double u = vel_x[Scalar::index(x, y)];
            double v = vel_y[Scalar::index(x, y)];

            double r = rho[Scalar::index(x, y)];

            // Initialise equilibrium
            double vel_sqr = u * u + v * v;

            for (int q = 0; q < Q; q++) {
                double vel_dot = u * c[q][0] + v * c[q][1];
                f1[Field::index(x, y, q)] = weights[q] * r * (
                    1.
                    + 3.0 * vel_dot
                    - 1.5 * vel_sqr
                    + 4.5 * vel_dot * vel_dot
                );
            }
        }
    }

    delete vel_x;
    delete vel_y;
    delete rho;

    for (int t = 0; t < val; t++) {
        // Collision step
        for (int y = 0; y < size_y; y++) {
            for (int x = 0; x < size_x; x++) {
                // -- Calculate macroscopic properties --

                double rho = 0.0;  // Density
                double u = 0.0;    // Velocity x
                double v = 0.0;    // Velocity y

                for (int q = 0; q < Q; q++) {
                    double dist = f1[Field::index(x, y, q)];

                    rho += dist;          // 0th order moment
                    u += c[q][0] * dist;  // 1st order moment
                    v += c[q][1] * dist;  // ^^^
                }
                // Velocity = Momentum / Density
                u /= rho;
                v /= rho;

                // Initialise equilibrium
                double vel_sqr = u * u + v * v;

                for (int q = 0; q < Q; q++) {
                    double vel_dot = u * c[q][0] + v * c[q][1];
                    double feq = weights[q] * rho * (
                        1.
                        + 3.0 * vel_dot
                        - 1.5 * vel_sqr
                        + 4.5 * vel_dot * vel_dot
                    );
                    f1[Field::index(x, y, q)] = (1 - omega) * f1[Field::index(x, y, q)] + omega * feq;
                }
            }
        }

        // Streaming step
        for (int y = 0; y < size_y; y++) {
            for (int x = 0; x < size_x; x++) {
                // Boundary types do not have distributions
                if (type_lattice[Scalar::index(x, y)] != PointType::Fluid)
                    continue;

                // Standard streaming w/ periodic boundary
                for (int q = 0; q < Q; q++) {
                    int xn = x + c[q][0];
                    int yn = y + c[q][1];

                    // Wrap coordinates over domain
                    if (xn >= size_x | xn < 0)
                        xn = mod(xn, size_x);
                    if (yn >= size_y | yn < 0)
                        yn = mod(yn, size_y);

                    switch (type_lattice[Scalar::index(xn, yn)]) {
                        case PointType::Fluid: {
                            f2[Field::index(xn, yn, q)] = f1[Field::index(x, y, q)];
                            break;
                        }
                        case PointType::NoSlip: {
                            f2[Field::index(x, y, bounce_back[q])] = f1[Field::index(x, y, q)];
                            break;
                        }
                        case PointType::Slip_V: {
                            f2[Field::index(x, yn, bounce_forward_v[q])] = f1[Field::index(x, y, q)];
                            break;
                        }
                        case PointType::Slip_H: {
                            f2[Field::index(xn, y, bounce_forward_h[q])] = f1[Field::index(x, y, q)];
                            break;
                        }
                    }
                }
            }
        }

        // Swap lattices to repeat process
        double* temp = f1;
        f1 = f2;
        f2 = temp;
    }

    mat2hdf5(f1, size_x, size_y);

    delete f1;
    delete f2;
    delete type_lattice;

    return 0;
}