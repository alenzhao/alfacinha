# Copyright 2007 Laurent GUEGUEN
# 
# This file is part of alfacinha..
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

DIR = alfacinha
SRCDIR = SrcCPP
MODDIR = Modules

#######################
#   cibles
#######################

.PHONY: all clean clean_all src exec tgz

all:
	$(MAKE) -C Modules all

clean:
	-rm $(MODDIR)/*.so
	-rm $(MODDIR)/*.pyc
	-cd $(SRCDIR); rm *.o; rm -r Pythondev

clean_all:
	-rm $(MODDIR)/*.so
	-rm $(MODDIR)/[^_]*.py
	-rm $(MODDIR)/*.pyc
	-cd $(SRCDIR); rm *.o; rm -r Pythondev; rm *.cxx

src:
	cd ..; tar -cvf $(DIR)/alfacinha_src.tar $(DIR)/COPYING $(DIR)/makefile \
	$(DIR)/*.py  $(DIR)/$(MODDIR)/*.py $(DIR)/$(MODDIR)/makefile \
	$(DIR)/$(SRCDIR)/*.inl $(DIR)/$(SRCDIR)/*.h $(DIR)/$(SRCDIR)/*.cpp \
	$(DIR)/$(SRCDIR)/*.cxx $(DIR)/$(SRCDIR)/makefile 

exec:
	cd ..; tar -cvf $(DIR)/alfacinha.tar $(DIR)/*.py $(DIR)/Modules/*.py \
	$(DIR)/$(MODDIR)/*.so

exgz:
	cd ..; tar -cvf $(DIR)/alfacinha.tar $(DIR)/*.py $(DIR)/Modules/*.py \
	$(DIR)/$(MODDIR)/*.so; gzip -f $(DIR)/alfacinha.tar

tgz:
	cd ..; tar -cvf $(DIR)/alfacinha_src.tar $(DIR)/COPYING $(DIR)/makefile \
	$(DIR)/*.py  $(DIR)/$(MODDIR)/*.py $(DIR)/$(MODDIR)/makefile \
	$(DIR)/$(SRCDIR)/*.inl $(DIR)/$(SRCDIR)/*.h $(DIR)/$(SRCDIR)/*.cpp \
	$(DIR)/$(SRCDIR)/*.cxx $(DIR)/$(SRCDIR)/makefile; \
	gzip -f $(DIR)/alfacinha_src.tar

ogz:
	cd ..; tar -cvf $(DIR)/alfacinha_orig.tar $(DIR)/COPYING $(DIR)/makefile \
	$(DIR)/*.py  $(DIR)/$(MODDIR)/*.py $(DIR)/$(MODDIR)/makefile \
	$(DIR)/$(SRCDIR)/*.inl $(DIR)/$(SRCDIR)/*.h $(DIR)/$(SRCDIR)/*.cpp \
	$(DIR)/$(SRCDIR)/*.i $(DIR)/$(SRCDIR)/makefile; \
	gzip -f $(DIR)/alfacinha_orig.tar

