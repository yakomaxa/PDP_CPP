#include "ShortSegmentRemover.hpp"
#include "PDPParameters.hpp"

void ShortSegmentRemover::cleanup(std::vector<Domain>& domains) {
    int ndom = (int)domains.size();

    for (int j = 0; j < ndom; j++) {
        int n = 0;
        bool allshort = true;
        // count the length of segments for this domain.
        for (int i = 0; i < domains[j].getNseg(); i++) {
            int seglen = domains[j].getSegmentAtPos(i).getTo() - domains[j].getSegmentAtPos(i).getFrom() + 1;
            if (seglen >= 30) allshort = false;
            n += seglen;
        }

        if (n < PDPParameters::MIN_DOMAIN_LENGTH || allshort) {
            ndom--;
            domains.erase(domains.begin() + j);
            j--;
        }
    }
}

