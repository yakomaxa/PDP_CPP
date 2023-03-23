#ifndef CUT_DOMAIN_HPP
#define CUT_DOMAIN_HPP

#include <vector>
#include <iostream>
#include "Atom.hpp"
#include "PDPDistanceMatrix.hpp"
#include "Domain.hpp"
#include "CutSites.hpp"
#include "CutValues.hpp"
#include "Cut.hpp"

class CutDomain {
    private:
        int ndom;
        std::vector<Domain> domains;
        std::vector<std::vector<int>> dist;
        std::vector<Atom> ca;

    public:
       static bool verbose;
       CutDomain(std::vector<Atom> ca, PDPDistanceMatrix pdpMatrix);
       void cutDomain(Domain dom,
                      CutSites cut_sites,
                      PDPDistanceMatrix pdpMatrix);
       std::vector<Domain> getDomains();
    
};

#endif  // CUT_DOMAIN_HPP

