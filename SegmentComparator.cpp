#include "SegmentComparator.hpp"

int SegmentComparator::operator()(const Segment& v1, const Segment& v2) const {
    //return v1.compareTo(v2);
    int result = 0;
    if (v1 > v2){
        result = 1;
    };
    if (v1 < v2){
        result = -1;
    };
    if (v1 == v2){
        result = 0;
    };
    return result;
}

