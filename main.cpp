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

int main(int argc, char *argv[]){
  std::string filename=argv[1];
  Structure s = Structure(filename);
  std::vector<Domain> domains;

    std::vector<Atom> ca = s.getRepresentativeAtomArray();

    GetDistanceMatrix distMaxCalculator;
    printf("A\n");
    PDPDistanceMatrix pdpMatrix = distMaxCalculator.getDistanceMatrix(ca);
    printf("B\n");
    Domain dom;
    //Chain c = ca[0].getGroup().getChain();
    //dom.setId("D"+c.getStructure().getPDBCode()+c.getId()+"1");
    dom.setId("testDomain");
    dom.setSize((int)ca.size());
    dom.setNseg(1);
    dom.getSegmentAtPos(0).setFrom(0);
    dom.getSegmentAtPos(0).setTo(int(ca.size())-1);
    printf("C\n");
    CutSites* cutSites = new CutSites();
    
    printf("%i\n",cutSites->cut_sites.size());
    // Do the initial splitting
    CutDomain cutDomain(ca,pdpMatrix);
    printf("D\n");
    cutDomain.cutDomain(dom, *cutSites, pdpMatrix);

    printf("E\n");
    domains =  cutDomain.getDomains();
    printf("F\n");
    listdomains(domains);
    printf("G\n");
    // Cluster domains
    domains = ClusterDomains::cluster(domains, pdpMatrix);
    printf("H\n");
    listdomains(domains);
    // Remove short segments
    ShortSegmentRemover::cleanup(domains);
    printf("FINAL!!\n");
    listdomains(domains);

    return 0;
}

