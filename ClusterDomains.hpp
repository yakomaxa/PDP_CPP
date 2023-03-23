#ifndef CLUSTER_DOMAINS_HPP
#define CLUSTER_DOMAINS_HPP

#include <vector>
#include "Domain.hpp"
#include "PDPDistanceMatrix.hpp"

class ClusterDomains {
public:
    ClusterDomains();
    static long getTotalContacts(std::vector<Domain> &domains,PDPDistanceMatrix& pdpDistMatrix, Domain& d1,Domain& d2);
    static std::vector<Domain> cluster(std::vector<Domain>& domains,PDPDistanceMatrix& pdpDistMatrix);
private:
    static int ndom;
    
};

#endif // CLUSTER_DOMAINS_HPP

