#ifndef PDP_DISTANCE_MATRIX_H
#define PDP_DISTANCE_MATRIX_H

#include <vector>
#include "Atom.hpp"
#include "PDPDistanceMatrix.hpp"

class GetDistanceMatrix{
public:
    PDPDistanceMatrix getDistanceMatrix(std::vector<Atom> protein);
};

#endif

