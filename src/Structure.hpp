//
//  Structure.hpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//

#ifndef Structure_hpp
#define Structure_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <gemmi/mmread.hpp>
#include <zlib.h>
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/util.hpp>

#include <stdio.h>
#include <vector>
#include <string>
#include <zlib.h>
#include <iomanip>

#include <iostream>
#include <fstream>

#include "Atom.hpp"

class Structure{
private:

public:
  gemmi::Structure structure;
  int numResidues;
  Structure(std::string filename);
  std::vector<Atom> getRepresentativeAtomArray();
  std::vector<int> tailofchain;
  void writePDB();
};

#endif /* Structure_hpp */
