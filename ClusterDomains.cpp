#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "ClusterDomains.hpp"
#include "PDPDistanceMatrix.hpp"
#include "Domain.hpp"
#include "CutDomain.hpp"
#include "PDPParameters.hpp"

static bool verbose = CutDomain::verbose;
static int ndom;

ClusterDomains::ClusterDomains(){
    
}

long calc_S(const int a1,
            const int b1,
            const int a2,
            const int b2,
            PDPDistanceMatrix pdpDistMatrix) {

    long contacts = 0;
    std::vector<std::vector<int>> dist = pdpDistMatrix.getDist();

    for (int i = a1; i <= b1; i++) {
        for (int j = a2; j <= b2; j++) {
            contacts += dist[i][j];
        }
    }

    return contacts;
};

long  ClusterDomains::getTotalContacts(std::vector<Domain>& domains,PDPDistanceMatrix& pdpDistMatrix,Domain& i,Domain& j) {
            long total_contacts = 0;
        for(int k=0; k<i.getNseg(); k++) {
            for(int l=0; l<j.getNseg(); l++) {
                long contacts = calc_S( j.getSegmentAtPos(l).getFrom(), j.getSegmentAtPos(l).getTo(), i.getSegmentAtPos(k).getFrom(), i.getSegmentAtPos(k).getTo(), pdpDistMatrix);
                total_contacts += contacts;
            }
        }
        return total_contacts;
};
    
std::vector<Domain> combine(std::vector<Domain> &domains, int Si, int Sj, double maximum_value) {
    if (verbose)
        std::cout << "  +++  combining domains " << Si << " " << Sj << std::endl;
    
    std::vector<Domain> newdoms;
    
    for (int i = 0; i < domains[Sj].getNseg(); i++) {
        domains[Si].getSegmentAtPos(domains[Si].getNseg()).setFrom(domains[Sj].getSegmentAtPos(i).getFrom());
        domains[Si].getSegmentAtPos(domains[Si].getNseg()).setTo(domains[Sj].getSegmentAtPos(i).getTo());
        //domains[Si].nseg++;
        domains[Si].addNseg(1);
    }
    domains[Si].addSize(domains[Sj].getSize());
    
    int ndom = (int)domains.size();
    for (int i = 0; i < domains[ndom - 1].getNseg(); i++) {
        domains[Sj].getSegmentAtPos(i).setFrom(domains[ndom - 1].getSegmentAtPos(i).getFrom());
        domains[Sj].getSegmentAtPos(i).setTo(domains[ndom - 1].getSegmentAtPos(i).getTo());
    }
    
    for (int i = 0; i < domains.size(); i++) {
        if (i == Sj)
            continue;
        newdoms.push_back(domains[i]);
    }
    
    domains[Sj].setSize(domains[ndom - 1].getSize());
    domains[Sj].setNseg(domains[ndom - 1].getNseg());
    printf("hogehoge\n");
    ndom--;
    return newdoms;
};


std::vector<Domain> ClusterDomains::cluster(
                            std::vector<Domain>& domains,
                            PDPDistanceMatrix& pdpDistMatrix){
    
  printf("LINE000001\n");
    int ndom = (int)domains.size()-1;
  printf("LINE000002\n");
    int Si = -1;
    int Sj = -1;
    int Sis = -1;
    int Sjs = -1;
    int Sim = -1;
    int Sjm = -1;
    
    long total_max_contacts = 0;
    
    double maximum_values = PDPParameters::CUT_OFF_VALUE1S;
    double maximum_valuem = PDPParameters::CUT_OFF_VALUE1M;
    double maximum_value  = PDPParameters::CUT_OFF_VALUE1;
    
    if (ndom < 2) return domains;
    
    printf("NDOM=%i\n",ndom);
    do {
        for(int i=0;i<ndom-1;i++) {
            for(int j=i+1;j<ndom;j++) {
	      printf("i j =%i %i\n",i,j);
	      	printf("LIne00000XXXXXXXX\n");
                Domain d1 = domains[i];
	      	printf("LIne00000XXXXXXX2\n");
                Domain d2 = domains[j];
	      	printf("LIne00000XXXXXXX3\n");
                long total_contacts =  ClusterDomains::getTotalContacts(domains,pdpDistMatrix,d1,d2);
                std::cout << " pos: d1:" << i << " vs d2:" << j << " d1:" << d1.getSegmentAtPos(0).getFrom() << "-" << d1.getSegmentAtPos(0).getTo() << " " <<  d2.getSegmentAtPos(0).getFrom() << "-" << d2.getSegmentAtPos(0).getTo() << " " << total_contacts << std::endl;
                int size1dom1=domains[i].getSize();
                int size2dom2=domains[j].getSize();
                double minDomSize=std::min(size1dom1,size2dom2);
                double maxDomSize=std::max(size1dom1,size2dom2);
                printf("LIne00003\n");
                // set some limits on how big the domains can get
                if(minDomSize>150&&maxDomSize>1.5*minDomSize){
		  maxDomSize=1.5*minDomSize;
		}else if(maxDomSize>2*minDomSize){
		  maxDomSize=2*minDomSize;
		}
		printf("LIne00004\n");
                long size1= std::min(PDPParameters::MAXSIZE,(int)minDomSize);
                long size2= std::min(PDPParameters::MAXSIZE,(int)maxDomSize);
                minDomSize=std::min(pow(minDomSize,1.6/3)+PDPParameters::RG1,pow(minDomSize,1.4/3)+pow(PDPParameters::TD1,1.6/3)+PDPParameters::RG1);
                maxDomSize=std::min(pow(maxDomSize,1.6/3)+PDPParameters::RG1,pow(maxDomSize,1.4/3)+pow(PDPParameters::TD1,1.6/3)+PDPParameters::RG1);
                
		printf("LIne00005\n");
                total_max_contacts=(long)(minDomSize*maxDomSize*10);
                if(size1>130){
		  total_max_contacts=(long)(minDomSize*maxDomSize*9);
		}
		printf("LIne00006\n");
                double S_value= total_contacts/(double)total_max_contacts;
                
                if (S_value  > maximum_value) {
                    maximum_value = S_value;
                    Si = i;
                    Sj = j;
                };
		printf("LIne00007\n");
                
                if (S_value > maximum_valuem&&size1<70) {
                    maximum_valuem = S_value;
                    Sim = i;
                    Sjm = j;
                };
		printf("LIne00008\n");
                if (S_value > maximum_values&&size1<52) {
                    maximum_values = S_value;
                    Sis = i;
                    Sjs = j;
                };
		printf("LIne00009\n");
                total_contacts = 0;
                
                total_max_contacts = 0;
                
            }
	    printf("LIne00010\n");            
        }
        

        if (maximum_value > PDPParameters::CUT_OFF_VALUE1) {
            domains = combine(domains,Si, Sj, maximum_value);
            maximum_value = PDPParameters::CUT_OFF_VALUE1-.1;
            maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
            maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
            
        }else if (maximum_valuem > PDPParameters::CUT_OFF_VALUE1M) {
            domains = combine(domains, Sim, Sjm, maximum_valuem);
            maximum_value =  PDPParameters::CUT_OFF_VALUE1-.1;
            maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
            maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
        }else if (maximum_values > PDPParameters::CUT_OFF_VALUE1S) {
            
            domains = combine(domains, Sis, Sjs, maximum_values);
            maximum_value = PDPParameters::CUT_OFF_VALUE1-.1;
            maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
            maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
            
        }else {
            
            maximum_value = -1.0;
            maximum_values = -1.0;
            maximum_valuem = -1.0;
            
        }
        
    } while (maximum_value > 0.0 || maximum_values > 0.0 || maximum_valuem > 0.0);
    return domains;
};


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

