#ifndef INDEX_H
#define INDEX_H

class Field {
   private:
    static int shape_x;
    static int shape_y;
    static int shape_q;

   public:
    static int size;
    static void set_field_shape(int shape_x, int shape_y, int shape_q);
    static unsigned int index(int x, int y, int q);
};

class Scalar {
   private:

   public:
    static int shape_x;
    static int shape_y;
    static int size;
    static void set_scalar_shape(int shape_x, int shape_y);
    static unsigned int index(int x, int y);
};

#endif