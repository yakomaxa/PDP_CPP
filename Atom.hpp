#ifndef ATOM_HPP
#define ATOM_HPP

#include <stdio.h>
#include <string>

class Atom {
private:
    float x, y, z;
    
public:
    Atom();
    ~Atom();
    void setX(float x);
    void setY(float y);
    void setZ(float z);
    float getX() ;
    float getY() ;
    float getZ() ;
};

#endif

