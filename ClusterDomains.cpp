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
int ClusterDomains::ndom;
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
        domains[Si].addNseg(1);
    }   
    domains[Si].addSize(domains[Sj].getSize());
    
    //int ndom = (int)domains.size();
    for (int i = 0; i < domains[ClusterDomains::ndom - 1].getNseg(); i++) {
      domains[Sj].getSegmentAtPos(i).setFrom(domains[ClusterDomains::ndom - 1].getSegmentAtPos(i).getFrom());
      domains[Sj].getSegmentAtPos(i).setTo(domains[ClusterDomains::ndom - 1].getSegmentAtPos(i).getTo());
    }
    
    for (int i = 0; i < domains.size(); i++) {
        if (i == Sj)
            continue;
        newdoms.push_back(domains[i]);
    }
    
    domains[Sj].setSize(domains[ClusterDomains::ndom - 1].getSize());
    domains[Sj].setNseg(domains[ClusterDomains::ndom - 1].getNseg());
    //printf("hogehoge\n");
    ClusterDomains::ndom--;
    return newdoms;
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




std::vector<Domain> ClusterDomains::cluster(
                            std::vector<Domain>& domains,
                            PDPDistanceMatrix& pdpDistMatrix){
    

    ClusterDomains::ndom = (int)domains.size();
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


    if(verbose){
      listdomains(domains);
    }
    
    if (ClusterDomains::ndom < 2){
      return domains;
    }

    
    do {
        for(int i=0;i<ClusterDomains::ndom-1;i++) {
            for(int j=i+1;j<ClusterDomains::ndom;j++) {
                Domain d1 = domains.at(i);
                Domain d2 = domains.at(j);
                long total_contacts =  ClusterDomains::getTotalContacts(domains,pdpDistMatrix,d1,d2);
                std::cout << " pos: d1:" << i << " vs d2:" << j << " d1:" << d1.getSegmentAtPos(0).getFrom() << "-" << d1.getSegmentAtPos(0).getTo() << " " <<  d2.getSegmentAtPos(0).getFrom() << "-" << d2.getSegmentAtPos(0).getTo() << " " << total_contacts << std::endl;
                int size1dom1=domains[i].getSize();
                int size2dom2=domains[j].getSize();
                double minDomSize=std::min(size1dom1,size2dom2);
                double maxDomSize=std::max(size1dom1,size2dom2);
                // set some limits on how big the domains can get
                if(minDomSize>150&&maxDomSize>1.5*minDomSize){
		  maxDomSize=1.5*minDomSize;
		}else if(maxDomSize>2*minDomSize){
		  maxDomSize=2*minDomSize;
		}
                long size1= std::min(PDPParameters::MAXSIZE,(int)minDomSize);
                long size2= std::min(PDPParameters::MAXSIZE,(int)maxDomSize);
                minDomSize=std::min(pow(minDomSize,1.6/3)+PDPParameters::RG1,pow(minDomSize,1.4/3)+pow(PDPParameters::TD1,1.6/3)+PDPParameters::RG1);
                maxDomSize=std::min(pow(maxDomSize,1.6/3)+PDPParameters::RG1,pow(maxDomSize,1.4/3)+pow(PDPParameters::TD1,1.6/3)+PDPParameters::RG1);
                
                total_max_contacts=(long)(minDomSize*maxDomSize*10);
                if(size1>130){
		  total_max_contacts=(long)(minDomSize*maxDomSize*9);
		}
                double S_value= total_contacts/(double)total_max_contacts;
		if(verbose) {
		  printf(" size1=%d size2=%d minDomSize=%5.2f maxDomSize=%5.2f total_contacts = %d \n", size1,size2,minDomSize,maxDomSize,total_contacts);
		  printf(" total_contacts = %d total_max_contacts = %d\n", total_contacts, total_max_contacts);
		  printf(" maximum_value = %f S_value = %f\n",maximum_value, S_value);
		}
                if (S_value  > maximum_value) {
                    maximum_value = S_value;
                    Si = i;
                    Sj = j;
                };
                
                if (S_value > maximum_valuem&&size1<70) {
                    maximum_valuem = S_value;
                    Sim = i;
                    Sjm = j;
                };

                if (S_value > maximum_values&&size1<52) {
                    maximum_values = S_value;
                    Sis = i;
                    Sjs = j;
                };

                total_contacts = 0;
                
                total_max_contacts = 0;
                
            }
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

