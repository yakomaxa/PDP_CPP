#include "CutDomain.hpp"

bool CutDomain::verbose = true;
int ndom;

CutDomain::CutDomain(std::vector<Atom> ca, PDPDistanceMatrix pdpMatrix){
    dist = pdpMatrix.getDist();
    this->ca = ca;
    ndom = 0;
    domains = std::vector<Domain>();
}

void CutDomain::cutDomain(Domain dom, CutSites cut_sites, PDPDistanceMatrix pdpMatrix) {
 
    int i, site;
    
    Domain dom1 = Domain();
    Domain dom2 = Domain();
    
    CutValues val = CutValues();
    val.s_min = 100;
    val.site2 = 0;
    val.first_cut = true;
    
    Cut cut = Cut();
    
    site = cut.cut(ca, dom, val, dist, pdpMatrix);
    
    if (site < 0) {
        domains.push_back(dom);
        dom.setScore(val.s_min);
        ndom++;
        return;
    }
    
    
//cut_sites.cut_sites[cut_sites.ncuts++] = site;
//cut_sites.cut_sites[cut_sites.ncuts++] = site;
    cut_sites.addNcuts(1);
    //cut_sites.setNcuts(cut_sites.getNcuts(),site);
    cut_sites.cut_sites[cut_sites.getNcuts()]=site;
    
    dom1.setSize(0);
    dom1.setNseg(0);
    dom2.setSize(0);
    dom2.setNseg(0);
    if (val.site2 == 0) {
        for (i = 0; i < dom.getNseg(); i++) {
            if (site > dom.getSegmentAtPos(i).getTo()) {
                dom1.getSegmentAtPos(dom1.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
                dom1.getSegmentAtPos(dom1.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
                dom1.addNseg(1);
                dom1.addSize(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
            } else if (site < dom.getSegmentAtPos(i).getFrom()) {
                dom2.getSegmentAtPos(dom2.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
                dom2.getSegmentAtPos(dom2.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
                dom2.addNseg(1);
                dom2.addSize(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
            } else if (site == dom.getSegmentAtPos(i).getFrom()) {
                dom2.getSegmentAtPos(dom2.getNseg()).setTo(site-1);
                dom2.getSegmentAtPos(dom2.getNseg()).setFrom(site-1);
                dom2.addNseg(1);
                dom1.getSegmentAtPos(dom1.getNseg()).setTo(site);
                dom1.getSegmentAtPos(dom1.getNseg()).setFrom(site);
                dom1.addNseg(1);
                dom1.addSize(1);
                dom2.addSize(1);
            } else if (site == dom.getSegmentAtPos(i).getTo()) {
                dom1.getSegmentAtPos(dom1.getNseg()).setTo(site);
                dom1.getSegmentAtPos(dom1.getNseg()).setFrom(site);
                dom1.addNseg(1);
                dom2.getSegmentAtPos(dom2.getNseg()).setTo(site+1);
                dom2.getSegmentAtPos(dom2.getNseg()).setFrom(site+1);
                dom2.addNseg(1);
                dom1.addSize(1);
                dom2.addSize(1);
            } else if (site > dom.getSegmentAtPos(i).getFrom() && site < dom.getSegmentAtPos(i).getTo()) {
                dom1.getSegmentAtPos(dom1.getNseg()).setTo(site);
                dom1.getSegmentAtPos(dom1.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
                dom1.addNseg(1);
                dom2.getSegmentAtPos(dom2.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
                dom2.getSegmentAtPos(dom2.getNseg()).setFrom(site+1);
                dom2.addNseg(1);
                dom1.addSize(site - dom.getSegmentAtPos(i).getFrom() + 1);
                dom2.addSize(dom.getSegmentAtPos(i).getTo() - site + 1);
            }
        }
    }
};

std::vector<Domain> CutDomain::getDomains(){
    return domains;
};



