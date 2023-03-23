#include "Structure.hpp"
#include "Atom.hpp"
//#include "Chain.hpp"
#include "GetDistanceMatrix.hpp"
#include "PDPDistanceMatrix.hpp"
#include "CutSites.hpp"
#include "CutDomain.hpp"
#include "ClusterDomains.hpp"
#include "ShortSegmentRemover.hpp"
#include "Domain.hpp"
#include "PDPParameters.hpp"
#include "LocalDomainParser.hpp"

static void listdomains(std::vector<Domain>& domains) {
  int i = -1;
  for (Domain& dom : domains) {
    i++;
    std::cout << "DOMAIN:" << i << " size:" << dom.getSize() << " " << dom.getScore() << std::endl;
    std::vector<Segment>& segments = dom.getSegments();

    for (Segment& s : segments) {
      std::cout << "   Segment: " << s << std::endl;
    }
  }
};

int main(){
    
    Structure s = Structure("./5cw9.cif");
    std::vector<Domain> domains;

    std::vector<Atom> ca = s.getRepresentativeAtomArray();

    GetDistanceMatrix distMaxCalculator;
    PDPDistanceMatrix pdpMatrix = distMaxCalculator.getDistanceMatrix(ca);

    Domain dom;
    //Chain c = ca[0].getGroup().getChain();
    //dom.setId("D"+c.getStructure().getPDBCode()+c.getId()+"1");
    dom.setId("testDomain");
    dom.setSize((int)ca.size());
    dom.setNseg(1);
    dom.getSegmentAtPos(0).setFrom(0);
    dom.getSegmentAtPos(0).setTo(int(ca.size())-1);

    CutSites cutSites;

    // Do the initial splitting
    CutDomain cutDomain(ca,pdpMatrix);
    cutDomain.cutDomain(dom, cutSites, pdpMatrix);
    domains =  cutDomain.getDomains();
    listdomains(domains);
    // Cluster domains
    domains = ClusterDomains::cluster(domains, pdpMatrix);
    listdomains(domains);
    // Remove short segments
    ShortSegmentRemover::cleanup(domains);

    listdomains(domains);

    return 0;
}

