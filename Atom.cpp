//
//  Atom.cpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//

#include "Atom.hpp"

Atom::Atom(){
    
};

Atom::~Atom() {
    // Nothing to do here
}

void Atom::setX(float x){
    this->x = x;
};

void Atom::setY(float y){
    this->y = y;
};

void Atom::setZ(float z){
    this->z = z;
};


float Atom::getX(){
    return this->x;
} ;

float Atom::getY(){
    return this->y;

} ;
float Atom::getZ(){
    return this->z;
};

