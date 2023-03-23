#ifndef LocalDomainParser_hpp
#define LocalDomainParser_hpp

#include <stdio.h>
#include <vector>
#include "Structure.hpp"
#include "ClusterDomains.hpp"

class LocalDomainParser{
public:
    std::vector<Domain> suggestDomains(Structure s);
};

#endif /* LocalDomainParser_hpp */
