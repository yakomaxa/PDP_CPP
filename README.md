# PDP++
Under development.

Protein Domain Parser (PDP) is a structure-baased protein domain parser, which cuts the protein chain into domain by maximizing the intradomain contacts and minimizing the interdomain contacts. I love this software, so I ported Java implementation found in BioJava into C++.

Though the behavior is under validation, my implementation has become faster than orignal binary. It also can cut multi-chain protein structure into domain, and also read mmCIF in addition to PDB file (Thanks to gemmi).  

Original paper:
PDP: protein domain parser [https://academic.oup.com/bioinformatics/article/19/3/429/258369]

Java ported version this C++ implementation is based on:
https://github.com/biojava/biojava/tree/master/biojava-structure/src/main/java/org/biojava/nbio/structure/domain

Now considering give this LGPL license after Java-ported PDP source code.

This software uses gemmi for structure-reading. This part follows gemmi's license (The Mozilla Public License).

pdppymol.py, a pymol extention to call PDP, is under MIT-license.


TODO

import author information from original java code for appropriate parts
