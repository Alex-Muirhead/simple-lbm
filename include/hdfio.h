#include <map>
#include <string>

#include "hdf5.h"

using namespace std;

#ifndef DATA_IO_H
#define DATA_IO_H

class InputData {
   public:
    hid_t file;
    unsigned int shape_x;
    unsigned int shape_y;
    unsigned int timesteps;
    unsigned int savestep;

    InputData(hid_t file);
    static InputData open(const char *filename);
    void close();
    void read_scalar(const char *dataname, double* buffer);
    void configure_domain(int* buffer);
};

class OutputData {
   private:
    map<string, hid_t> datasets;

   public:
    hid_t file;

    static OutputData open(const char *filename);
    static OutputData create(const char *filename);
    void close();
    void write_scalar(const char *dataname, double *data);
    void append_scalar(const char *dataname, double *data);
};

#endif