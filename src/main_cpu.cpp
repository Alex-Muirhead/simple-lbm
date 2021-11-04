#include <array>
#include <iostream>
#include <sstream>

#include "constants.h"
#include "hdfio.h"
#include "index.h"

using namespace std;

enum PointType : int {
    Fluid = 1 << 0,   // Standard streaming
    NoSlip = 1 << 1,  // Bounce-back boundary
    Slip_V = 1 << 2,  // Vertical bounce-forward boundary
    Slip_H = 1 << 3,  // Horiztonal bounce-forward boundary
};

static unsigned int size_x;
static unsigned int size_y;

void init_eq(double* f, int* phase, double* rho, double* vel_x, double* vel_y) {
    for (unsigned int y = 0; y < size_y; y++) {
        for (unsigned int x = 0; x < size_x; x++) {
            if (phase[Scalar::index(x, y)] != PointType::Fluid) {
                // Set boundary lattice at zero
                for (int q = 0; q < Q; q++) {
                    f[Field::index(x, y, q)] = 1.0;  // Initial density
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
                f[Field::index(x, y, q)] = weights[q] * r * (1. + 3.0 * vel_dot - 1.5 * vel_sqr + 4.5 * vel_dot * vel_dot);
            }
        }
    }
}

void collide(double* f, double omega) {
    // Loop over the lattice dimensions
    for (unsigned int y = 0; y < size_y; y++) {
        for (unsigned int x = 0; x < size_x; x++) {
            // -- Calculate macroscopic properties --

            double r = 0.0;  // Density
            double u = 0.0;  // Velocity x
            double v = 0.0;  // Velocity y

            for (int q = 0; q < Q; q++) {
                double dist = f[Field::index(x, y, q)];

                r += dist;            // 0th order moment
                u += c[q][0] * dist;  // 1st order moment
                v += c[q][1] * dist;  // ^^^
            }
            // Velocity = Momentum / Density
            u /= r;
            v /= r;

            // Precompute velocity norm
            double vel_sqr = u * u + v * v;

            // Can unroll loop if necessary
            for (int q = 0; q < Q; q++) {
                double vel_dot = u * c[q][0] + v * c[q][1];
                double feq = weights[q] * r * (1. + 3.0 * vel_dot - 1.5 * vel_sqr + 4.5 * vel_dot * vel_dot);
                f[Field::index(x, y, q)] = (1 - omega) * f[Field::index(x, y, q)] + omega * feq;
            }
        }
    }
}

void stream(double* f_src, double* f_dst, int* phase) {
    for (unsigned int y = 0; y < size_y; y++) {
        for (unsigned int x = 0; x < size_x; x++) {
            // Boundary types do not have distributions
            if (phase[Scalar::index(x, y)] != PointType::Fluid)
                continue;

            // Standard streaming w/ periodic boundary
            for (unsigned int q = 0; q < Q; q++) {
                // Wrap coordinates over domain
                unsigned int xn = int(size_x + x + c[q][0]) % size_x;
                unsigned int yn = int(size_y + y + c[q][1]) % size_y;

                unsigned int next_index;

                switch (phase[Scalar::index(xn, yn)]) {
                    case PointType::NoSlip: {
                        next_index = Field::index(x, y, bounce_back[q]);
                        break;
                    }
                    case PointType::Slip_V: {
                        next_index = Field::index(x, yn, bounce_forward_v[q]);
                        break;
                    }
                    case PointType::Slip_H: {
                        next_index = Field::index(xn, y, bounce_forward_h[q]);
                        break;
                    }
                    default: {
                        next_index = Field::index(xn, yn, q);
                        break;
                    }
                }

                f_dst[next_index] = f_src[Field::index(x, y, q)];
            }
        }
    }
}

void calculate_flow_properties(double* f, double* rho, double* vel_x, double* vel_y) {
    for (unsigned int y = 0; y < size_y; y++) {
        for (unsigned int x = 0; x < size_x; x++) {
            // -- Calculate macroscopic properties --

            double r = 0.0;  // Density
            double u = 0.0;  // Velocity x
            double v = 0.0;  // Velocity y

            for (int q = 0; q < Q; q++) {
                double dist = f[Field::index(x, y, q)];

                r += dist;            // 0th order moment
                u += c[q][0] * dist;  // 1st order moment
                v += c[q][1] * dist;  // ^^^
            }
            // Velocity = Momentum / Density
            u /= r;
            v /= r;

            rho  [Scalar::index(x, y)] = r;
            vel_x[Scalar::index(x, y)] = u;
            vel_y[Scalar::index(x, y)] = v;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Expected filename" << endl;
        return -1;
    }

    double nu_lb = 0.092;                     // Lattice dynamic viscosity
    double omega = 1.0 / (3. * nu_lb + 0.5);  // Relaxation parameter

    InputData config = InputData::open(argv[1]);

    Field::set_field_shape(config.shape_x, config.shape_y, Q);
    Scalar::set_scalar_shape(config.shape_x, config.shape_y);

    size_x = config.shape_x;
    size_y = config.shape_y;

    int* type_lattice = new int[Scalar::size];

    config.configure_domain(type_lattice);

    double* rho   = new double[Scalar::size];
    double* vel_x = new double[Scalar::size];
    double* vel_y = new double[Scalar::size];

    config.read_scalar("initial/density", rho);
    config.read_scalar("initial/vel_x", vel_x);
    config.read_scalar("initial/vel_y", vel_y);

    double* f1 = new double[Field::size];
    double* f2 = new double[Field::size];

    // Set up equilibrium values for initial flow
    init_eq(f1, type_lattice, rho, vel_x, vel_y);

    OutputData output;
    output = OutputData::create("output.h5");
    output.write_scalar("density", rho);
    output.write_scalar("vel_x", vel_x);
    output.write_scalar("vel_y", vel_y);

    for (unsigned int t = 0; t < config.timesteps; t++) {
        // Collision step
        collide(f1, omega);

        // Streaming step
        stream(f1, f2, type_lattice);

        // Decide when to save / export the data
        if ((t + 1) % config.savestep == 0) {
            printf("Saving data from timestep %d\n", t + 1);
            calculate_flow_properties(f2, rho, vel_x, vel_y);
            output.append_scalar("density", rho);
            output.append_scalar("vel_x", vel_x);
            output.append_scalar("vel_y", vel_y);
        }

        // Swap lattices to repeat process
        double* temp = f1;
        f1 = f2;
        f2 = temp;
    }

    config.close();
    output.close();

    delete vel_x;
    delete vel_y;
    delete rho;

    delete f1;
    delete f2;
    delete type_lattice;

    return 0;
}