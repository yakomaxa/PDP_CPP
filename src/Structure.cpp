//
//  Structure.cpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//
#include "Structure.hpp"
#include "PDPParameters.hpp"

#include <gemmi/to_cif.hpp>
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_mmcif.hpp>
#define GEMMI_READ_CIF_IMPLEMENTATION
#include <gemmi/read_cif.hpp>

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
  int index=-1;
  int CA_flag=0;
  int chainid=0;
  int resi=0;
  int maxindex=0;
  for (gemmi::Model& model : this->structure.models){
    for (gemmi::Chain& chain : model.chains) {
      CA_flag=0;
      for (gemmi::Residue& residue : chain.residues) {
	for (gemmi::Atom &atom : residue.atoms) {
	  std::string elementname = atom.element.name();
	  if (atom.name == "CA" && elementname == "C"){
	    index += 1;
	    Atoms[index].setX(atom.pos.x);
	    Atoms[index].setY(atom.pos.y);
	    Atoms[index].setZ(atom.pos.z);
	    Atoms[index].setChain(chain.name);
	    resi=stoi(residue.seqid.str());
	    Atoms[index].setIndexOrg(resi);
	    Atoms[index].setChainId(chainid);
	    if (maxindex < resi){
	      maxindex = resi;
	    }	      
	    CA_flag=1;
	  }
	  if (atom.name == "CB" && elementname == "C"){
	    Atoms[index].setX(atom.pos.x);
	    Atoms[index].setY(atom.pos.y);
	    Atoms[index].setZ(atom.pos.z);
	  }
	}
      }
      if(CA_flag==1){
	this->tailofchain.push_back(index);
      }
      chainid++;
    }
  }
  PDPParameters::maxIndex = maxindex;
  return Atoms;
};


void Structure::writePDB(){
  gemmi::Structure st_tmp = this->structure;

  /*
  long long iatom = 0;
  gemmi::Structure st_tmp = st;
  for (gemmi::Model& model : st_tmp.models){
    for (gemmi::Chain& chain : model.chains) {
       for (gemmi::Residue& residue : chain.residues) {           
	 for (gemmi::Atom &atom : residue.atoms) {
	     iatom = iatom + 1;
	     atom.pos.x = (double)allatm[iatom].coord[1] ;
	     atom.pos.y = (double)allatm[iatom].coord[2] ;
	     atom.pos.z = (double)allatm[iatom].coord[3] ;	     
	 }
       }
    }
  }
  */
  
  /*std::string mmcifout = input.mmcifout;
  std::string imodel_str=std::to_string(imodel);
  if (input.dump_topN <= 0){
    mmcifout = mmcifout+"_MODEL"+imodel_str+".cif";
  }else{
    std::ostringstream ss;
    ss << std::setw(3) << std::setfill('0') << isub_out+1;
    std::string isub_out_str(ss.str());
    mmcifout = mmcifout+"_ranked"+isub_out_str+"_MODEL"+imodel_str+".cif";
  }
  std::ofstream os1(mmcifout);
  */
  std::string mmcifout = "testoutput.cif";
  std::ofstream os1(mmcifout);
  //generate cif document from Structure object
  gemmi::cif::Document doc = gemmi::make_mmcif_document(st_tmp);
  //object group
  gemmi::MmcifOutputGroups groups(true);
  groups.atoms = true;
  //use group_pdb i.e. ATOM and HETATM
  groups.group_pdb = true;
  // update the mmcif by Structure object
  gemmi::update_mmcif_block(st_tmp, doc.blocks[0], groups);
  //std::string file_q_str(input.file_q);
  gemmi::cif::write_cif_to_stream(os1, doc);
  
}
