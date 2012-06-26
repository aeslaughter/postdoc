#!/bin/bash
cd /home/slaughter/Documents/programs/build
make matlab-doc
cd ..
doxygen Doxyfile.web
cd /home/slaughter/Documents/web/postdoc
cp --recursive /home/slaughter/Documents/programs/doc/web/* ./doc/
