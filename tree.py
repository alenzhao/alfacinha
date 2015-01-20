# -*- coding: utf-8 -*-

"""Newick format parser module.

Manipulates trees in the Newick format. This format is mainly used to
describe phylogenetic binary trees but can have a much wider use. Some
details on the format can be found here:

http://evolution.genetics.washington.edu/phylip/newicktree.html

Examples of Newick format:

* (Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;

* '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'
        
Note that bracked delimited nested and unnested comments are ignored.

Sequences can be assigned to each node for simulating evolution.

Copyright 2007, Leonor Palmeira, <palmeira@biomserv.univ-lyon1.fr>.
"""

__author__ = "Leonor Palmeira <palmeira@biomserv.univ-lyon1.fr>"
__date__ = "17 October 2007"
__credits__ = """Guido van Rossum, for an excellent programming language."""

import string

import evol

#######################################################################
#######################################################################
########  Class Node

class Node:
    """Node is the smallest unit of definition for a tree.

    A Node is linked to its father-Node and its children-Nodes: it
    defines a sub-tree for which it is the root. All connected Nodes
    (including leaves and root) make up for the whole tree.
    """

    def __init__(self, **kw):
	"""Create a Node.
	
	Keyword argument to build directly the Node from a Newick
	string: [newick=string].

        Examples of Newick format:

        * '(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268,Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;'

        * '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'

        Note that bracked delimited nested and unnested comments are ignored.
    	"""

	self.__l=0 #length of the branch to the father-Node
	self.__lab="" #Node label
	self.__seq=None #Node sequence
	self.__boot=None #bootstrap value at the Node
	self.__father=None
	self.__children=[]
	if kw.has_key('newick'):
	    self.parser(kw['newick'])
	if kw.has_key('fic'):
	    self.read_nf(kw['fic'])


    def sequence(self):
        """Return the Sequence of the Node."""
        return self.__seq

    def lg(self):
        """Return the length of the edge to the father."""
        return self.__l

    def set_lg(self,l):
        """Set the length of the edge to the father to $2 if it is >=0."""
        if l>=0:
            self.__l=l

    def __getitem__(self,n):
        """Return the Node with this name."""
        for s in self.__children:
            if s.label()==n:
                return s
            else:
                t=s[n]
                if t:
                    return t
        return None
    
    def __str__(self):
	"""Return printable string in Newick format of the sub-tree
	defined by the Node.

"""
	return self.newick()

    def nb_children(self):
        """Return the number of children."""
        return len(self.__children)

    def label(self):
        """Return the label."""
        return self.__lab
    
    def read_nf(self,a):
	"""Read the $1 file containing a unique tree in Newick format
        and builds the Node from it.
"""

	if type(a)==str:
            f=open(a,'r')
            l=f.readline()
            f.close()
	    self.parser(l)
	else:
	    raise ValueError, "Invalid file name."
	return
	


    def write_newick(self,a, **kw):
	"""Write the Node in Newick format in $1 file.

	Keyword argument to determine whether 'append' or 'write'
	(overwrites) mode is turned on: [mode=string].
"""
	if kw.has_key('mode'):
	    mode=kw['mode']
	else:
	    raise ValueError, 'Should I append or overwrite file?'	
	if type(a)==str:
	    if mode=='append':
		f=open(a,'a')
	    elif mode=='write':
		f=open(a,'w')
            f.write(self.newick())
            f.close()
	else:
	    raise ValueError, "Unvalid file name. Unable to write."
	return
	
    def write_phylip(self,a, **kw):
	"""Write the Node in Phylip format in $1 file.

	Keyword argument to determine whether 'append' or 'write'
	(overwrites) mode is turned on: [mode=string].
"""
	if kw.has_key('mode'):
	    mode=kw['mode']
	else:
	    raise ValueError, 'Should I append or overwrite file?'	
	if type(a)==str:
	    if mode=='append':
		f=open(a,'a')
            elif mode=='write':
                f=open(a,'w')
	    f.write(self.phylip())
	    f.close()
	else:
	    raise ValueError, "Unvalid file name. Unable to write."
	return

    def write_matrix(self, a, **kw):
	"""Write the Node in simple format in $1 file.
	
	Keyword argument to determine whether 'append' or 'write'
	(overwrites) mode is turned on: [mode=string].
"""
	if kw.has_key('mode'):
	    mode=kw['mode']
	else:
	    raise ValueError, 'Should I append or overwrite file?'	
	s=self.matrix()+'\n'
	if type(a)==str:
	    if mode=='append':
		f=open(a,'a')
	    else:# mode=='write':
		f=open(a,'w')
	    f.write(s)
	    f.close()
	else:
	    raise ValueError, "Unvalid file name. Unable to write."
	return

#####################################
############## screen output methods:

    def newick(self):
	"""Return Newick string of the sub-tree defined by the Node.

	Follows the Newick standard for coding trees.
	Newick format example:
	       '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'
"""

	s=''
	end=''
	if self.__father==None:
	    end=';'
	i=0
	while i<len(self.__children):
	    if i==0:
		s+='('
	    s+=self.__children[i].newick()
	    if i!=(len(self.__children)-1):
		s+=','
	    else:
		s+=')'
	    i+=1
	if self.__boot!=None:
	    s+=str(self.__boot)
	elif self.__lab!=None:
	    s+=self.__lab
	if self.__l!=None:
	    s+=':'+str(self.__l)
	s+=end
	return s

    def phylip(self):
	"""Return output of sequences of the leaves defined by the
	Node in Phylip interleaved format.
	
"""
	
	s=''
	l=self.get_leaf_sequences()
	s+=str(len(l))+'\t'+str(len(l[0]))+'\n'
	a=[]
	for i in l:
	    this=''
	    name=string.join(i.name().split(),'_')
	    length=len(name)
	    if length>10:
		this+=i.name()[:10]
	    else:
		this+=i.name()+' '*(10-length)
	    seq=i.seq()
	    n=len(seq)/10
	    for j in range(n):
		this+=seq[j*10:j*10+10]+' '
	    this+=seq[n*10:]
	    a+=[this]
	for i in a:
	    s+=i[0:31]+'\n'
	m=len(a[0][31:])/21 # all lines are of the same length
	j=11
	while j<=m*21:
	    s+='\n'
	    for i in a:
		s+='          '+i[j+21:j+21*2]+'\n' # exactly 10 white spaces needed...
	    j+=22
	s+='\n'
	for i in a:
	    s+='          '+i[j+21:]+'\n'
	return s

    def matrix(self):
	"""Return additive distances matrix of the leaves.
	
"""

	d=self._matrix()
	n=d.keys()
	n.sort() #order lines and columns
	s='' #s=str(len(n))+"\n" # changed this from s=''
	for i in n:
#            s+=i+" " # changed this (added this line)
	    for j in n:
		s+=str(d[i][j])+'\t'
	    s+='\n'
	s=string.strip(s)
	return s
	
#####################################
############## miscellaneous methods:
	
    def get_leaf_labels(self):
	"""Return the list of labels of the leaves defined by the Node.

"""
	
	a=[]
	if self.__children!=[]:
	    for i in self.__children:
		a+=i.get_leaf_labels()
	else:
	    a+=[self.__lab]
	return a

    def get_leaf_sequences(self):
	"""Return the list of sequences of the leaves defined by the
	Node.

"""

	a=[]
	if self.__children!=[]:
	    for i in self.__children:
		a+=i.get_leaf_sequences()
	else:
	    a+=[self.__seq]
	return a

    def get_leaves(self):
	"""Return the list of leaves defined by the Node."""

	a=[]
	if self.__children!=[]:
	    for i in self.__children:
		a+=i.get_leaves()
	else:
	    a+=[self]
	return a

    def __imul__(self, x):
        """ Recursively multiplies the length of the branches
that are under this node by factor x (>0).
"""
        if x>0:
            for i in self.__children:
                i.__l*=x
                i*=x
        return self
        
    def __idiv__(self, x):
        """ Recursively divides the length of the branches
that are under this node by factor x (>0).
"""
        if x>0:
            for i in self.__children:
                i.__l/=x
                i/=x
        return self

#####################################
########### distance measure methods:

    def depth(self):
	"""Return the depth of this Node as its 'distance' (in number of branches) to the root.

WARNING: Branch lengths are not use.
"""
	
	d=0
	if self.__father:
	    d=self.__father.depth()+1
	return d

    def distance_root(self):
        """Return the distance (sum of branch lengths) from this Node
to the root.
"""

        if self.__father!=None:
            return self.__l+self.__father.distance_root()
        else:
            return self.__l
        
    def distance(self,n):
	"""Return the distance (sum of branch lengths) separating this
Node (leaf) from Node $1 (leaf).

"""

	if self==n:
	    d=float(0)
	else:
	    if self.depth()>n.depth():
		d=self.go_father().distance(n)+float(self.__l)
	    elif self.depth()<n.depth():
		d=n.go_father().distance(self)+float(n.__l)
	    elif self.depth()==n.depth():
		if self.__father!=n.__father:
		    d=self.__l+n.__l+self.go_father().distance(n.go_father())
		else:
		    d=self.__l+n.__l
	return d

    def _matrix(self):
	"""Return the distance matrix of all leaves, starting at the
root.
"""	
	self=self.go_root()
	l=self.get_leaves()
	d={}
	for i in range(len(l)):
	    d[l[i].__lab]={}
	    for j in range(len(l)):
		d[l[i].__lab][l[j].__lab]=l[i].distance(l[j])
	return d	
    
#####################################		
### 'walking along the tree' methods:

    def go_root(self):
	"""Return Root Node."""

	if self.__father:
	    self=self.__father.go_root()
	return self

    def go_father(self):
	"""Return father Node of this Node."""
	if self.__father:
	    self=self.__father
	else:
	    raise IndexError, 'Already at root. No more father to go to.'
	return self

    def go_child(self,n):
	"""Return $1 (number) child of this Node."""
	
	if self.__children!=[] and n<len(self.__children):
	    self=self.__children[n]
	else:
	    raise IndexError, 'Already at leaf. No more children to go to.'
	return self

#####################################
####### EvolSequence related methods:

    def assign_seq(self,seq):
	"""Assign a sequence to the Node.

	WARNING: Intended to be used with the EvolSequence class to
	simulate evolution.
"""
	
	self.__seq=seq
	return

    def evolve_seq(self,seq,mod, **kw): #verifier que les arguments marchent!!!
	"""Assign a sequence $1 to the Node, and evolve it along the
	sub-tree defined by the Node, according to model $2.

        WARNING: Must be used with the EvolSequence class in order to
        be able to simulate evolution.

	Keyword argument to determine how the number of substitutions
	on the sequence is approximated: [approx=string] 'Rounded' is
	the default. Keyword argument to determine the algorithm used
	to simulate neighbor-dependent substitutions: [algo=string]
	'Berard' is the default.

        Optional argument 'positions' is a list of all allowed
        positions. This argument defaults to all positions of the
        EvolSequence.
        
        WARNING: When an empty list is given, all positions are
        considered allowed.
"""

        self.__seq=seq.copy()
        self.__seq.g_name(self.__lab)
	self.__seq.evolve(mod,self.__l, **kw)
        
	for c in self.__children:
            c.evolve_seq(self.__seq,mod, **kw)
            
    def evolve_seg_seq(self,seq,dmod,**kw):
	"""Assign a sequence $1 to the Node, and evolve it along the
	sub-tree defined by the Node, according to models in
	dictionary $2.

        The dictionary items are (deb,fin):mod where the modele mod is
        applied to the positions in range [deb:fin]. If ranges
        overlap, a random model is applied on each substitution that
        occurs in the overlap.
        
        WARNING: Must be used with the EvolSequence class in order to
        be able to simulate evolution.

	Keyword argument to determine how the number of substitutions
	on the sequence is approximated: [approx=string] 'Rounded' is
	the default. Keyword argument to determine the algorithm used
	to simulate neighbor-dependent substitutions: [algo=string]
	'Berard' is the default.

        """

        self.__seq=seq.copy()
        self.__seq.g_name(self.__lab)
	self.__seq.evolve_seg(dmod,self.__l, **kw)
        
	for c in self.__children:
            c.evolve_seg_seq(self.__seq,dmod, **kw)
            
#####################################
#################### parsing methods:
    
    def clean(self,s):
	"""Clean the string in Newick format.

	Bracked delimited nested and unnested comments are eliminated.
"""
	s=s.strip()
	if s[len(s)-1]==';':
	    if s.count('(')!=s.count(')'):
		raise ValueError, "Opening parenthesis do not match closing parenthesis."
	    else:
		brackopen=s.count('[')
		brackclose=s.count(']')
		if brackopen!=0 or brackclose!=0:
		    if brackopen==brackclose:
			open=s.find('[')
			x=1
			for i in range(open+1,len(s)):
			    if x>0:
				if s[i]=='[':
				    x+=1
				elif s[i]==']':
				    x-=1
			    else:
				break
			s=s[0:open]+s[i:len(s)]
			s=self.clean(s)
		    else:
			raise ValueError, "Opening brackets do not match closing brackets."
	else:
	    raise ValueError, "Missing ';' at end."
	return s

    def _parser(self,s):
	"""Should not be directly used. Use parser() instead."""
	parenth=s.rfind(')')
	semicol=s.rfind(':')
	if semicol>parenth: # deal with node values
	    try:
		self.__l=float(s[semicol+1:len(s)])
	    except ValueError, e:
		raise ValueError, "Incorrect branch value -> must be numerical."
	    try:
		self.__boot=float(s[parenth+1:semicol])
	    except ValueError, e:
		self.__lab=s[parenth+1:semicol]
	elif semicol==-1 and parenth!=-1:
	    self.__lab=s[parenth+1:]
	elif semicol==-1:
	    raise ValueError, "Incorrect syntax."
	s=s[1:parenth] #eliminate external parenthesis
	x=0
	cuts=[0] #cutting points for nodes of the same depth
	for i in range(len(s)):
	    if s[i]=='(':
		x+=1
	    elif s[i]==')':
		x-=1
	    elif x==0 and s[i]==',':
		cuts+=[i+1]
	if cuts!=[0]:
	    cuts+=[len(s)]
	    i=0
	    while i<len(cuts)-1:
		child=s[cuts[i]:cuts[i+1]].strip(',')
		self.__children+=[Node()]
		for j in self.__children:
		    j.__father=self
		self.__children[len(self.__children)-1]._parser(child)
#		print "***********Going to the children*************"
		i+=1
	return
	
    def parser(self,s):
	"""Fill the Node's attributes from parsing $1 string.

	Follows the Newick standard for coding trees.
	Bracked delimited nested and unnested comments are ignored.
        
	Newick format example:
	       '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'

	WARNING: all interior nodes (if given) should be strings ->
	numerical values will be interpreted as the node's bootstrap
	value.
	"""

	s=self.clean(s)
        s="".join(s.split()).strip(";") # remplace s=s.strip(';')
	if s!=-1:
	    self._parser(s)
	else:
	    raise ValueError, "This should not have happened! Check behind your back."
        
