#ifndef SHORT_SEGMENT_REMOVER_HPP
#define SHORT_SEGMENT_REMOVER_HPP

#include <vector>
#include "Domain.hpp"

class ShortSegmentRemover {
public:
    static void cleanup(std::vector<Domain>& domains);
};

#endif /* SHORT_SEGMENT_REMOVER_HPP */

