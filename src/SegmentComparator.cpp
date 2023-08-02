#include "SegmentComparator.hpp"

bool SegmentComparator::operator()(const Segment& v1, const Segment& v2) const {
    //return v1.compareTo(v2);
  int result = false;
  //return 0;
  if (v1 > v2){
    result = false;
  };
  if (v1 < v2){
    result = true;
  };
  return result;
}

