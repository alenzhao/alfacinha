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

FSH = -shared -fPIC 
FSHD= -bundle -framework Python

OPS= -g

SRCDIR = ../SrcCPP

VPATH = .:$(SRCDIR)

######################

define faitmodule
@if [ `uname` = Darwin ]; then \
cd $(SRCDIR); set -x; $(CCC) $(OPS) $(FSHD) $(^F) -o ../Modules/$@ ; \
else \
cd $(SRCDIR); set -x; $(CCC) $(OPS) $(FSH)  $(^F) -o ../Modules/$@; \
fi
endef

############################
# MACROS
#############################

ELO=  descripteur.o #Liste.o
SEQUENCEO = sequence.o Fsequence.o random.o

#######################
#   cibles
#######################

all: _Csequencemodule.so _Fsequencemodule.so 

.PHONY: all

.PRECIOUS: %_wrap.cxx

%_wrap.cxx: %.i
	$(MAKE) -C $(SRCDIR) $(@F)

%_wrap.o: %_wrap.cxx 
	$(MAKE) -C $(SRCDIR) $(@F)

%.o:	%.cpp
	$(MAKE) -C $(SRCDIR) $(@F)

_Csequencemodule.so: $(addprefix $(SRCDIR)/, $(SEQUENCEO) sequence_wrap.o)
	$(faitmodule)

_Fsequencemodule.so: $(addprefix $(SRCDIR)/, Fsequence.o Fsequence_wrap.o random.o)
	$(faitmodule)


