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



int  PDPDistanceMatrix::getNclose_raw() {
    return nclose_raw;
}

void PDPDistanceMatrix::setNclose_raw(int nclose_raw) {
    this->nclose_raw = nclose_raw;
}

std::vector<int> PDPDistanceMatrix::getIclose_raw() {
    return iclose_raw;
}

void PDPDistanceMatrix::setIclose_raw(std::vector<int> iclose_raw) {
    this->iclose_raw = iclose_raw;
}

std::vector<int> PDPDistanceMatrix::getJclose_raw() {
    return jclose_raw;
}

void PDPDistanceMatrix::setJclose_raw(std::vector<int> jclose_raw) {
    this->jclose_raw = jclose_raw;
}

