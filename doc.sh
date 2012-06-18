#!/bin/bash
cd /home/slaughter/Documents/programs
doxygen Doxyfile.web
cd /home/slaughter/Documents/web/postdoc
cp --recursive /home/slaughter/Documents/programs/doc/web/html/* ./doc/
