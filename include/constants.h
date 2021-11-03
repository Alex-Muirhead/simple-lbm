static const int D = 2; // 2 dimensions to the simulation
static const int Q = 9; // 9 different velocity components in the distribution

static const float c[Q][D] = {
    // Stationary (v* = 0)
    {  0,  0 },
    // Cardinal directions (v* = 1)
    {  0, -1 },
    {  1,  0 },
    {  0,  1 },
    { -1,  0 },
    // Diagonal directions (v* = sqrt{2})
    {  1, -1 },
    {  1,  1 },
    { -1,  1 },
    { -1, -1 }
};

static const float weights[Q] = {
    // Stationary
    4./9,
    // Cardinal directions
    1./9,
    1./9,
    1./9,
    1./9,
    // Diagonal directions
    1./36,
    1./36,
    1./36,
    1./36
};

static const int bounce_back[Q] = {
    // Stationary
    /* 0 -> */ 0,
    // Cardinal directions
    /* 1 -> */ 3,
    /* 2 -> */ 4,
    /* 3 -> */ 1,
    /* 4 -> */ 2,
    // Diagonal directions
    /* 5 -> */ 7,
    /* 6 -> */ 8,
    /* 7 -> */ 5,
    /* 8 -> */ 6,
};

static const int bounce_forward_h[Q] = {
    // Stationary
    /* 0 -> */ 0,
    // Cardinal directions
    /* 1 -> */ 3,
    /* 2 -> */ 2,
    /* 3 -> */ 1,
    /* 4 -> */ 4,
    // Diagonal directions
    /* 5 -> */ 6,
    /* 6 -> */ 5,
    /* 7 -> */ 8,
    /* 8 -> */ 7,
};

static const int bounce_forward_v[Q] = {
    // Stationary
    /* 0 -> */ 0,
    // Cardinal directions
    /* 1 -> */ 1,
    /* 2 -> */ 4,
    /* 3 -> */ 3,
    /* 4 -> */ 2,
    // Diagonal directions
    /* 5 -> */ 8,
    /* 6 -> */ 7,
    /* 7 -> */ 6,
    /* 8 -> */ 5,
};