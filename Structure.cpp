//
//  Structure.cpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//

#include "Structure.hpp"

gemmi::Structure openCoordinatefile(const std::string & infilename){
  gemmi::MaybeGzipped inputfile(infilename);
  gemmi::CoorFormat infileformat = gemmi::coor_format_from_ext(inputfile.basepath());
  if(infileformat != gemmi::CoorFormat::Unknown){
    return gemmi::read_structure(inputfile, infileformat);
  }else{
    return gemmi::read_structure(inputfile, gemmi::CoorFormat::Pdb);
  }
};

Structure::Structure(std::string filename){
    this->structure=openCoordinatefile(filename);
    this->numResidues = 0;
    for (gemmi::Model& model : this->structure.models){
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                for (gemmi::Atom &atom : residue.atoms) {
                    std::string elementname = atom.element.name();
                    if (atom.name == "CA" && elementname == "C"){
                        this->numResidues += 1;
                    }
                }
            }
        }
    }
};

std::vector<Atom> Structure::getRepresentativeAtomArray(){
    std::vector<Atom> Atoms(this->numResidues);
    int index=0;
    for (gemmi::Model& model : this->structure.models){
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                for (gemmi::Atom &atom : residue.atoms) {
                    std::string elementname = atom.element.name();
		    if (atom.name == "CA" && elementname == "C"){
		      Atoms[index].setX(atom.pos.x);
		      Atoms[index].setY(atom.pos.y);
		      Atoms[index].setZ(atom.pos.z);
		      printf("CA\n");
		      index += 1;
                    }
		    if (atom.name == "CB" && elementname == "C"){
		      Atoms[index].setX(atom.pos.x);
		      Atoms[index].setY(atom.pos.y);
		      Atoms[index].setZ(atom.pos.z);
		      printf("CB\n");
                    }
                }
            }
        }
    }

    return Atoms;
};


