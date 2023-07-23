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
      total_contacts +=  contacts;
    }
  }
  return total_contacts;
};

bool  ClusterDomains::isContacting(Domain& i,Domain& j,std::vector<int> &iclose,std::vector<int> &jclose, int nclose){
  for(int k=0; k<i.getNseg(); k++) {
    int fromi=i.getSegmentAtPos(k).getFrom();
    int toi=i.getSegmentAtPos(k).getTo();
    for(int l=0; l<j.getNseg(); l++) {
      int fromj=j.getSegmentAtPos(l).getFrom();
      int toj=j.getSegmentAtPos(l).getTo();
      for (int n = 0 ; n < nclose ; n++){
	if (fromi <= iclose[n] && toi >= iclose[n] &&
	    fromj <= jclose[n] && toj >= jclose[n] 
	    ){
	  return true;
	}
	if (fromi <= jclose[n] && toi >= jclose[n] &&
	    fromj <= iclose[n] && toj >= iclose[n] 
	    ){
	  return true;
	}
      }                  
    }
  }
  printf("NOCONTACT\n");
  return false;
};
    
std::vector<Domain> combine(std::vector<Domain> &domains, int Si, int Sj, double maximum_value) {
  if (verbose){
        std::cout << "  +++  combining domains " << Si << " " << Sj << std::endl;
  }
  printf("Hoge01\n")  ;
  std::vector<Domain> newdoms;
  printf("Hoge02\n")  ;
  for (int i = 0; i < domains[Sj].getNseg(); i++) {
    domains[Si].getSegmentAtPos(domains[Si].getNseg()).setFrom(domains[Sj].getSegmentAtPos(i).getFrom());
    domains[Si].getSegmentAtPos(domains[Si].getNseg()).setTo(domains[Sj].getSegmentAtPos(i).getTo());    
    domains[Si].addNseg(1);  
  }
  printf("Hoge02_2\n")  ;
  for (int k : domains[Sj].getContacted()){
    printf("Contact %i\n",k);
    domains[Si].removeContacted(k);
    domains[Si].pushbackContacted(k);
  }
  printf("Hoge03\n")  ;
  domains[Si].addSize(domains[Sj].getSize());
  printf("Hoge04\n")  ;
  ///int ndom = (int)domains.size();
  for (int i = 0; i < domains[ClusterDomains::ndom - 1].getNseg(); i++) {  
    domains[Sj].getSegmentAtPos(i).setFrom(domains[ClusterDomains::ndom - 1].getSegmentAtPos(i).getFrom());
    domains[Sj].getSegmentAtPos(i).setTo(domains[ClusterDomains::ndom - 1].getSegmentAtPos(i).getTo());
  }
  printf("Hoge05\n")  ;


  
  for (int i = 0; i < ClusterDomains::ndom; i++) {
    printf("Hoge05_1\n")  ;
    if (i == Sj){
      continue;
    }
    printf("Hoge05_2\n");
    newdoms.push_back(domains[i]);
    printf("Hoge05_3\n")  ;
  }

  printf("Hoge06\n")  ;
  domains[Sj].setSize(domains[ClusterDomains::ndom - 1].getSize());
  domains[Sj].setNseg(domains[ClusterDomains::ndom - 1].getNseg());
  printf("hogehoge\n");
  ClusterDomains::ndom--;
  return newdoms;
};


std::vector<Domain> ClusterDomains::cluster(
                            std::vector<Domain>& domains,
                            PDPDistanceMatrix& pdpDistMatrix){
    
  printf("LINE000001\n");
  ClusterDomains::ndom = (int)domains.size();
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
  
  if (ClusterDomains::ndom < 2){
    return domains;
  }
  
  Domain d1;
  Domain d2;
  std::vector<std::vector<int>> total_contacts_list(ClusterDomains::ndom, std::vector<int>(ClusterDomains::ndom));
  
  std::vector<int> i_can_contact(ClusterDomains::ndom*ClusterDomains::ndom);
  std::vector<int> j_can_contact(ClusterDomains::ndom*ClusterDomains::ndom);
  int n_can=0;
  
  std::vector<int> icontacted(ClusterDomains::ndom*ClusterDomains::ndom);
  std::vector<int> jcontacted(ClusterDomains::ndom*ClusterDomains::ndom);
  int n_contact_pair=0;
  int i,j;
  
  int nclose_raw = pdpDistMatrix.getNclose_raw();
  std::vector<int> iclose_raw = pdpDistMatrix.getIclose_raw();
  std::vector<int> jclose_raw = pdpDistMatrix.getJclose_raw();
  printf("%i\n",nclose_raw);
  for(int i=0;i<ClusterDomains::ndom-1;i++) {
    d1 = domains.at(i);
    for(int j=i+1;j<ClusterDomains::ndom;j++) {
      d2 = domains.at(j);
      if (ClusterDomains::isContacting(d1,d2,iclose_raw,jclose_raw,nclose_raw)){
	i_can_contact[n_can]=i;
	j_can_contact[n_can]=j;	    
	n_can++;
      }
    }
  }
  
  for(int n = 0 ; n < n_can; n++) {
    i = i_can_contact[n];
    j = j_can_contact[n];
    
    if (i < j){
      
    }else{
      int tmp = i;
      i = j;
      j = tmp;	
    }
    d1 = domains.at(i);    
    d2 = domains.at(j);
    long total_contacts=ClusterDomains::getTotalContacts(domains,pdpDistMatrix,domains[i],domains[j]);
    std::cout << " pos: d1:" << i << " vs d2:" << j << " d1:" << d1.getSegmentAtPos(0).getFrom() << "-" << d1.getSegmentAtPos(0).getTo() << " " <<  d2.getSegmentAtPos(0).getFrom() << "-" << d2.getSegmentAtPos(0).getTo() << " " << total_contacts_list[i][j] << std::endl;
    if (total_contacts > 0){
      domains[i].pushbackContacted(j);
      domains[j].pushbackContacted(i);
      n_contact_pair+=1;
    }
  }
  
  printf("CLUSTERDOMAIN::NDOM=%i\n",ClusterDomains::ndom);
  do {
    printf("DO\n");
    printf("CLUSTERDOMAIN::NDOM=%i\n",ClusterDomains::ndom);
    //    std::vector<std::vector<int>> contacts;
    for(int i=0;i<ClusterDomains::ndom;i++) {
      for (int j : domains[i].getContacted()){
	if (j==i){
	  continue;
	}
	long total_contacts = ClusterDomains::getTotalContacts(domains,pdpDistMatrix,domains[i],domains[j]);
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
	
	total_max_contacts=(long)(minDomSize*maxDomSize*10);
	if(size1>130){
	  total_max_contacts=(long)(minDomSize*maxDomSize*9);
	}
	printf("LIne00006\n");
	
	double S_value= total_contacts/(double)total_max_contacts;
	printf("S_value %f\n,", S_value);
	
	
      if(verbose) {
	printf("size1=%i size2=%i minDomSize=%f maxDomSize=%f total_contacts = %i \n", size1,size2,minDomSize,maxDomSize,total_contacts);
	printf(" total_contacts = %i total_max_contacts = %i\n", total_contacts, total_max_contacts);
	printf(" maximum_value = %f S_value = %f\n",maximum_value, S_value);
      }
      
      if (S_value  > maximum_value) {
	maximum_value = S_value;
	Si = i;
	Sj = j;
      };
      printf("LIne00007\n");
      
      if (S_value > maximum_valuem && size1<70) {
	maximum_valuem = S_value;
	Sim = i;
	Sjm = j;
      };
      printf("LIne00008\n");
      if (S_value > maximum_values && size1<52) {
	maximum_values = S_value;
	Sis = i;
	Sjs = j;
      };
      printf("LIne00009\n");
      total_contacts = 0;
      
      total_max_contacts = 0;
      
      printf("LIne00010\n");            
      }
    }

    printf("LIne00011\n");            
    if (maximum_value > PDPParameters::CUT_OFF_VALUE1) {
      printf("LIne00012\n");            
      domains = combine(domains,Si, Sj, maximum_value);
      printf("LIne00013\n");            
      maximum_value = PDPParameters::CUT_OFF_VALUE1-.1;
      maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
      maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
      
    }else if (maximum_valuem > PDPParameters::CUT_OFF_VALUE1M) {
      printf("LIne00014\n");            
      domains = combine(domains, Sim, Sjm, maximum_valuem);
      printf("LIne00015\n");            
      maximum_value =  PDPParameters::CUT_OFF_VALUE1-.1;
      maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
      maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
    }else if (maximum_values > PDPParameters::CUT_OFF_VALUE1S) {
      printf("LIne00016\n");            
      domains = combine(domains, Sis, Sjs, maximum_values);
      printf("LIne00017\n");            
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

