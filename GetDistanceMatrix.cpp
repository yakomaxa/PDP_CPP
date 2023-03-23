#include <iostream>
#include <math.h>
#include "Group.hpp"
#include "Atom.hpp"
#include "PDPDistanceMatrix.hpp"
#include "PDPParameters.hpp"
#include "GetDistanceMatrix.hpp"

float getDistance(Atom atom1,Atom atom2);
float getDistance(Atom atom1,Atom atom2){
    float dist;
    float dx;
    float dy;
    float dz;
    
    dx = atom1.getX() - atom2.getX();
    dy = atom1.getY() - atom2.getY();
    dz = atom1.getZ() - atom2.getZ();
    
    dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    return dist;
}


PDPDistanceMatrix GetDistanceMatrix::getDistanceMatrix(std::vector<Atom> protein)
{
    
    
    std::vector<std::vector<int>> dist(PDPParameters::MAXLEN+3,
                                  std::vector<int>(PDPParameters::MAXLEN+3));
    int i,j;
    float d,dt1,dt2,dt3,dt4;
    int nclose=0;
    std::vector<int> iclose(PDPParameters::MAXLEN*PDPParameters::MAXLEN);
    std::vector<int> jclose(PDPParameters::MAXLEN*PDPParameters::MAXLEN);
    
    if(sizeof(protein)/sizeof(protein[0]) >= PDPParameters::MAXLEN) {
        std::cerr << sizeof(protein)/sizeof(protein[0]) << " protein.len > MAXLEN " << PDPParameters::MAXLEN << std::endl;
        //return 0;
        exit;
    }
    for(i=0; i < sizeof(protein)/sizeof(protein[0]); i++) {
        for(j=i; j < sizeof(protein)/sizeof(protein[0]); j++) {
            dist[i][j]=0;
            dist[j][i]=0;
            d=0;
            
            Atom ca1 = protein.at(i);
            Atom ca2 = protein.at(j);
            //Group g1 = ca1.getGroup();
            //Group g2 = ca2.getGroup();
            
            //Atom* cb1 = getCBeta(g1);
            //Atom* cb2 = getCBeta(g2);
            //Atom cb1 = NULL;
            //Atom cb2 = NULL;
            //bool hasCbeta1 = cb1 != NULL;
            //bool hasCbeta2 = cb2 != NULL;
            bool hasCbeta1=false;
            bool hasCbeta2=false;
            
            dt1=81;
            dt2=64;
            dt3=49;
            dt4=36;
            
            if(hasCbeta1 && hasCbeta2) {
//                float distance = Calc::getDistance(cb1,cb2);
//                d+= distance*distance;
            }
            else if(hasCbeta1 && !hasCbeta2) {
//                float distance = 999;
//                distance = Calc::getDistance(cb1, ca2);
//                d += distance * distance;
            }
            else if(!hasCbeta1&&hasCbeta2) {
//                float distance = Calc::getDistance(ca1, cb2);
//                d += distance * distance;
            }
            else if(!hasCbeta1&&!hasCbeta2) {
                float distance = getDistance(ca1, ca2);
                d += distance * distance;
            }
            
            if(d<dt1) {
                dist[i][j]=1;
                dist[j][i]=1;
                if(d<dt2) {
                    dist[i][j]=2;
                    dist[j][i]=2;
                    if(j-i>35) {
                        iclose[nclose]=i;
                        jclose[nclose]=j;
                        nclose++;
                    }
                    if(d<dt3) {
                        dist[i][j]=4;
                        dist[j][i]=4;
                        if(d<dt4) {
                            dist[i][j]=6;
                            dist[j][i]=6;
                        }
                    }
                }
            }
        }
    }
    
    for(int i=1;i<protein.size();i++) {
        for(int j=i;j<protein.size()-1;j++) {
            if(dist[i][j]>=2&&j-i>5) {
                if((dist[i-1][j-1]>=2&&dist[i+1][j+1]>=2)||(dist[i-1][j+1]>=2&&dist [i+1][j-1]>=2))  {
                    dist[i][j]+=4;
                    dist[j][i]+=4;
                }
                else if(i>2&&j<protein.size()-2) {
                    if((dist[i-3][j-3]>=1&&dist[i+3][j+3]>=1)||(dist[i-3][j+3]>=1&&dist[i+3][j-3]>=1)) {
                        dist[i][j]+=4;
                        dist[j][i]+=4;
                    }
                    else if(i>3&&j<protein.size()-3) {
                        if(((dist[i-3][j-3]>=1||dist[i-3][j-4]>=1||dist[i-4][j-3]>=1||dist[i-4][j-4]>=1)&&
                            (dist[i+4][j+4]>=1||dist[i+4][j+3]>=1||dist[i+3][j+3]>=1||dist[i+3] [j+4]>=1))
                           ||((dist[i-4][j+4]>=1||dist[i-4][j+3]>=1||dist[i-3][j+4]>=1||dist[i-3][j+3]>=1)&&
                              (dist[i+4][j-4]>=1||dist[i+4][j-3]>=1||dist[i+3][j-4]>=1||dist[i+3][j-3]>=1))) {
                            dist[i][j]+=4;
                            dist[j][i]+=4;
                            
                        }
                    }
                }
            }
        }
    }
    PDPDistanceMatrix matrix;
    matrix.setNclose(nclose);
    std::vector iv = iclose;
    matrix.setIclose(iv);
    matrix.setJclose(jclose);
    matrix.setDist(dist);
    return matrix;
}


