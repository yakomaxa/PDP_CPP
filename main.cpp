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
//#include "LocalDomainParser.hpp"

static void listdomains(std::vector<Domain>& domains) {
  int i = -1;
  for (Domain& dom : domains) {
    i++;
    std::cout << "create DOMAIN" << i << ", resi ";
    std::vector<Segment>& segments = dom.getSegments();

    int flag=0;
    for (Segment& s : segments) {
      if (flag>0){
	std::cout << "+" ;
      }
      std::cout << s ;
      flag++;
    }
    std::cout << ";" << std::endl;
  }
};

int main(int argc, char *argv[]){
  std::string filename=argv[1];
  printf("---------Reading structure\n");
  Structure s = Structure(filename);
  printf("---------Reading structure Done\n");
  std::vector<Domain> domains;

  printf("---------Get Repr atoms\n");
  std::vector<Atom> ca = s.getRepresentativeAtomArray();
  printf("---------Get Repr atoms Done\n");
  GetDistanceMatrix distMaxCalculator;

  printf("---------distMat creation\n");
  PDPDistanceMatrix pdpMatrix = distMaxCalculator.getDistanceMatrix(ca);
  printf("---------distMat creation Done\n");

  Domain dom;
  //Chain c = ca[0].getGroup().getChain();
  //dom.setId("D"+c.getStructure().getPDBCode()+c.getId()+"1");
  printf("---------Setting domain info\n");
  dom.setId("testDomain");
  dom.setSize((int)ca.size());
  dom.setNseg(1);

  dom.getSegmentAtPos(0).setFrom(0);
  dom.getSegmentAtPos(0).setTo(int(ca.size())-1);
  
  CutSites cutSites = CutSites();
  printf("---------Setting domain info done\n");  
  // Do the initial splitting
  printf("---------Initial splitting\n");
  std::vector<int> init_cutsites;
  // add the head residue of the chains except for the first one
  s.tailofchain.erase(s.tailofchain.end()-1);
  for (int i : s.tailofchain){
    printf("ADDING INITAIAL SITE%i\n",i);
    init_cutsites.push_back(i+1);
  }
  CutDomain cutDomain(ca,pdpMatrix, init_cutsites);
  cutDomain.cutDomain(dom, cutSites, pdpMatrix, init_cutsites);
  printf("---------Initial splitting done\n");  


  domains =  cutDomain.getDomains();
  listdomains(domains);
  // Cluster domains
  printf("---------Clustering domains\n");  
  domains = ClusterDomains::cluster(domains, pdpMatrix);
  printf("---------Clustering domains Done\n");  
  listdomains(domains);
  // Remove short segments
  printf("---------Cleanup \n");  
  ShortSegmentRemover::cleanup(domains);
  printf("---------Cleanup Done\n");  
  printf("FINAL!!\n");
  listdomains(domains);
  
    return 0;
}

