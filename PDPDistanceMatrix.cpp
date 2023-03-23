#include "PDPDistanceMatrix.hpp"

PDPDistanceMatrix::PDPDistanceMatrix() {
    // Default constructor
}

std::vector<std::vector<int>> PDPDistanceMatrix::getDist() {
    return dist;
}

void PDPDistanceMatrix::setDist(std::vector<std::vector<int>> dist) {
    this->dist = dist;
}

int  PDPDistanceMatrix::getNclose() {
    return nclose;
}

void PDPDistanceMatrix::setNclose(int nclose) {
    this->nclose = nclose;
}

std::vector<int> PDPDistanceMatrix::getIclose() {
    return iclose;
}

void PDPDistanceMatrix::setIclose(std::vector<int> iclose) {
    this->iclose = iclose;
}

std::vector<int> PDPDistanceMatrix::getJclose() {
    return jclose;
}

void PDPDistanceMatrix::setJclose(std::vector<int> jclose) {
    this->jclose = jclose;
}

