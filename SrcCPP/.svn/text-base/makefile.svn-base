# Copyright 2004 Laurent GUEGUEN
# 
# This file is part of Sarment.
# 
# Sarment is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# Sarment is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Sarment; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#############################
## C++

CCC= g++

OPS = -g

INC = $(shell \
	if [ `uname` = Darwin ]; then if [ -x /Library/Frameworks/Python.framework/Headers ] ; then echo '/Library/Frameworks/Python.framework/Headers'; elif [-x /System/Library/Frameworks/Python.framework/Headers] ; then echo '/System/Library/Frameworks/Python.framework/Headers'; fi \
	elif [ -x /usr/include/python2.5 ] ; then echo '/usr/include/python2.5'; \
	elif [ -x /usr/include/python2.4 ] ; then echo '/usr/include/python2.4'; \
	elif [ -x /usr/include/python2.3 ] ; then echo '/usr/include/python2.3'; \
	elif [ -x /usr/include/python2.2 ]; then echo '/usr/include/python2.2'; \
	elif [ -x /appli2/Python/include/python2.2 ]; then echo '/appli2/Python/include/python2.2'; \
	elif [ -x /usr/local/include/python2.2 ]; then echo '/usr/local/include/python2.2'; \
	fi)

#Uncomment 
SPECIFIC = $(shell if [ `uname` = Darwin ]; then echo '-Dunix'; \
	else echo '-DPyMODINIT_FUNC=void'; fi)

FPY = -fPIC -IPythondev $(SPECIFIC)

FLAGS = -felide-constructors -Wunused -Winline $(FPY)

##########################"
## SWIG

SWIG = swig
SFLAGS= -c++ -python -shadow

MODDIR = ../Modules

VPATH = .:$(SRCDIR)

#############################
# cibles
############################

.PRECIOUS: %_wrap.cxx

%_wrap.cxx: %.i
	$(SWIG) $(SFLAGS) $^ 
	-mv *.py $(MODDIR)

%_wrap.o: %_wrap.cxx Pythondev
	if [ -x Pythondev ]; then  $(CCC) $(OPS) $(FPY) -c $(<F); fi

%.o: %.cpp 
	$(CCC) $(OPS) $(FLAGS) -c $(<F)

Pythondev:
ifneq ($(INC), )
	ln -s $(INC) Pythondev
else
	@echo 
	@echo "**********************************"
	@echo "ERROR: No directory link for SrcCPP/Pythondev"
	@echo "**********************************"
	@echo 
endif


