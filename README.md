# PDP++
Under development.

Protein Domain Parser (PDP) is a structure-based algorithm/program for splitting protein into domains, which cuts the protein chain into domains by maximizing the intradomain contacts and minimizing the interdomain contacts. I love this software, so I ported a Java implementation found in BioJava into C++.

Though the behavior is under validation, my implementation has become faster than orignal binary. It also can cut multi-chain protein structure into domains, and also read mmCIF in addition to PDB file (Thanks to gemmi).  

Original paper:
PDP: protein domain parser [https://academic.oup.com/bioinformatics/article/19/3/429/258369]

Java ported version this C++ implementation is based on:
https://github.com/biojava/biojava/tree/master/biojava-structure/src/main/java/org/biojava/nbio/structure/domain

Now considering give this LGPL license after Java-ported PDP source code.

This software uses gemmi for structure-reading. This part follows gemmi's license (The Mozilla Public License).

pdppymol.py, a pymol extention to call PDP, is under MIT-license.


TODO
(0) Priority: Let the binary to dump parsed structure files
(1) Import author information from original java code for appropriate parts
(2) Make systematic check routine
