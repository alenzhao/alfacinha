# -*- coding: utf-8 -*-
"""Evolutionary models with neighbor-dependent substitutions module.

Defines evolutionary models of evolution. This class inherits most of
its functions from Porportion class of compte module.

Copyright 2007, Leonor Palmeira, <palmeira@biomserv.univ-lyon1.fr>.
"""

__author__ = "Leonor Palmeira <palmeira@biomserv.univ-lyon1.fr>"
__date__ = "17 October 2007"
__credits__ = """Guido van Rossum, for an excellent programming language."""

import re

from compte import *

#######################################################################
#######################################################################
########  Class Model

class Model(Proportion):
    """A Model describes substitution processes acting on a sequence.

    A Model is any kind of model which can incorporate:
    
    * any simple substitutions desired
    * any neighbor-dependent substitutions (note that multiple substitutions are not implemented).

    It should be built from a file in the following format: (where A|A
    is not allowed, and CG|TG is an incorrect syntax).

    A|T   6
    A|G   5
    [... for simple substitution rates]

    Cg|T 10 [... for neighbor-dependent substitutions due to the right
    neighbor. Here C is substituted by T]
    cG|A 8 [... for neighbor-dependent substitutions due to the left
    neighbor. Here C is substituted by A]
    Cga|T 6 [... for neighbor-dependent substitutions due to the two
    right neighboring bases. Here C is substituted by T]

    Note that the neighboring bases affecting a given substitution can
    be as large as desired.
    """

    def __init__(self, **kw):
	"""Create a Model.
	
	Keyword argument to build directly the Model from a file in
	the appropriate format: [fic=string].
    	"""
	
	Proportion.__init__(self,**kw)#fic=kw['fic'])
	self.__a=self.alph_upper() #upper-case letters only
	self.__next_all=self.next_all()
	self.__max=self.max_subst()
        
    def next_extended(self,s):
	"""Return the dictionary of {prefix:[[postfix, value],
	[postfix, value]]} which has $1 in the prefix.

	$1 is a letter.
	"""
	
	list=[]
        for j in self.prefixes():
            k=j.find(s)
            if k!=-1:
                l=self.next(j)
                for a in l:
                    b=[j[:k].upper(),j[k+1:].upper(),a[0],a[1]]
                    list.append(b)
	return list

    def next_all(self):
	"""Return a dictionary of {alpha:{prefix:[postfix,value]}} of
	all possible capital letters and the possible substitutions
	associated.
	
	"""
	
	all={}
	a=self.__a
	for s in a:
	    all[s]=self.next_extended(s)
	return all
	
    def sum_subst(self,s):
	"""Return the sum of substitutions that can act on $1.

	$1 can be a sub-word.
	"""
	
	sum=0
	n=self.next_extended(s)
	for l in n:
            sum+=l[3]
	return sum
    
    def max_subst(self):
	"""Return maximum number of substitutions for all letters from
	the alphabet.
	"""

	s=0
	a=self.__a
        for i in a:
            j=self.sum_subst(i)
            if s<=j:
                s=j
	return s

    def alph_upper(self):
        """Return ONLY the capital letters contained in the alphabet.

        """

	a=[]
	for i in self.alph():
	    if i.isupper():
		a+=i
        return a
