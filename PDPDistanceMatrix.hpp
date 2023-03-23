#ifndef PDPCORE_PDPDISTANCEMATRIX_H
#define PDPCORE_PDPDISTANCEMATRIX_H

#include <vector>
#include <stdio.h>

class PDPDistanceMatrix {
private:
    std::vector<std::vector<int>> dist;
    int nclose;
    std::vector<int> iclose;
    std::vector<int> jclose;
    
public:
    PDPDistanceMatrix();
    std::vector<std::vector<int>> getDist();
    void setDist(std::vector<std::vector<int>> dist);
    
    int getNclose();
    void setNclose(int nclose);
    
    std::vector<int> getIclose();
    void setIclose(std::vector<int> iclose);
    std::vector<int> getJclose();
    void setJclose(std::vector<int> jclose);
};

#endif // PDPCORE_PDPDISTANCEMATRIX_H


