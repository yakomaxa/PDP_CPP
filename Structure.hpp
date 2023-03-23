//
//  Structure.hpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//

#ifndef Structure_hpp
#define Structure_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <zlib.h>

#include "Atom.hpp"
#include "./gemmi/mmread.hpp"
#include "./gemmi/gz.hpp"
#include "./gemmi/cif.hpp"
#include "./gemmi/util.hpp"

class Structure{
private:

public:
    gemmi::Structure structure;
    int numResidues;    
    Structure(std::string filename);
    std::vector<Atom> getRepresentativeAtomArray();
};

#endif /* Structure_hpp */
