#include "index.h"

void Field::set_field_shape(int shape_x, int shape_y, int shape_q) {
    Field::shape_x = shape_x;
    Field::shape_y = shape_y;
    Field::shape_q = shape_q;
    Field::size = shape_x * shape_y * shape_q;
}

unsigned int Field::index(int x, int y, int q) {
    return shape_q * (shape_x * y + x) + q;
}

int Field::shape_x;
int Field::shape_y;
int Field::shape_q;
int Field::size;


void Scalar::set_scalar_shape(int shape_x, int shape_y) {
    Scalar::shape_x = shape_x;
    Scalar::shape_y = shape_y;
    Scalar::size = shape_x * shape_y;
}

unsigned int Scalar::index(int x, int y) {
    return shape_x * y + x;
}

int Scalar::shape_x;
int Scalar::shape_y;
int Scalar::size;