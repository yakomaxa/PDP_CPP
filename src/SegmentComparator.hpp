#ifndef SEGMENT_COMPARATOR_HPP
#define SEGMENT_COMPARATOR_HPP

#include <iostream>
#include <vector>
#include "Segment.hpp"

class SegmentComparator {
public:
    bool operator()(const Segment& v1, const Segment& v2) const;
};

#endif /* SEGMENT_COMPARATOR_HPP */

