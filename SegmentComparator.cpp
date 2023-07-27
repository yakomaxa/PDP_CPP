#include "SegmentComparator.hpp"

int SegmentComparator::operator()(const Segment& v1, const Segment& v2) const {
    //return v1.compareTo(v2);
    int result = 0;
    //printf("WOOOOOOOOOOOOOOOOOOOo\n");
    if (v1 > v2){
        result = 1;
	//printf("COND01\n");
    };
    if (v1 < v2){
        result = -1;
	//printf("COND02\n");
    };
    if (v1 == v2){
        result = 0;
	//printf("COND03\n");
    };
    return result;
}

