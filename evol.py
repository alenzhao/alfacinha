# -*- coding: utf-8 -*-
"""Sequence evolution simulator module.

Defines sequences able to evolve. This class inherits most of
its functions from Sequence class of compte sequence.

Copyright 2007, Leonor Palmeira, <palmeira@biomserv.univ-lyon1.fr>.
"""

__author__ = "Leonor Palmeira <palmeira@biomserv.univ-lyon1.fr>"
__date__ = "17 October 2007"
__credits__ = """Guido van Rossum, for an excellent programming language."""

from compte import *
from sequence import *

#from misc import *

#######################################################################
#######################################################################
########  Class EvolSequence

class EvolSequence(Sequence):
    """An EvolSequence is a Sequence with the ability to evolve.

    It can be built from a file in the appropriate format. Only fasta
    format (with .fa or .fst extension) is allowed.
    """

#####################################
######### replacing Sequence methods:

#    def __init__(self, **kw): #n'est plus utile
#        sequence.Sequence.__init__(self, **kw)
#        ls=len(self)

    def copy(self):
        """Return a NEW copy EvolSequence of the EvolSequence.
	"""
	
        g=EvolSequence()
	g._Seq__gen.copie(self._Seq__c_elem())
	return g

#####################################
#####################################
##### EvolSequence specific methods:

#     def write_fasta(self,a, **kw):
# 	"""Write the EvolSequence in FASTA format in $1 file.

# 	Keyword argument to determine whether 'append' or 'write'
# 	(overwrites) mode is turned on: [mode=string].
# 	"""
	
# 	if kw.has_key('mode'):
# 	    mode=kw['mode']
# 	else:
# 	    raise ValueError, 'Should I append or overwrite file?'	
# 	if type(a)==str:
# 	    if mode=='append':
# 		f=open(a,'a')
# 	    elif mode=='write':
# 		f=open(a,'w')
#             else:
#                 raise ValueError, 'Should I append or overwrite file?'
                
#             f.write(self.fasta())
#             f.close()
# 	else:
# 	    raise ValueError, "Unvalid file name. Unable to write."
# 	return

    def compare(self,gen):
	"""Return the observed differences between this EvolSequence
	and sequence $1.
	"""
        ls=len(self)
	lg=len(gen)
	if ls==lg:
	    j=0
	    for i in range(ls):
		if self[i]!=gen[i]:
		    j+=1
	else:
	    raise ValueError, "Sequences of different lengths."
	return j

    def pick_position(self, segments=range(0)):
        """Return a random position in the EvolSequence.

        Optional argument 'segments' is a list of all allowed
        segments. Otherwise a segment is a SORTED list [beg,end], with
        end position excluded. The probability of a position to be
        chosen is proportional to the number of segments this position
        belongs to. This argument defaults to the whole sequence.
        """

        if segments==[]:
            l=len(self)
        else:
            l=0
            for s in segments:
                l+=s[1]-s[0]
        i=random.randint(0,l-1)
        if segments==[]:
            return i
        d=0
        while i>=(segments[d][1]-segments[d][0]):
           i-=segments[d][1]-segments[d][0]
           d+=1
        return i+segments[d][0]

           
#####################################
######## methods for some statistics:

    def freq(self): #gerer les ^ ?
	"""Return the observed frequencies of nucleo tides in the
	EvolSequence.
        
	"""
	
	c=Compte()
	c.add_seq(self,1)
	c.__idiv__(float(len(self)))
	return c
   
    def difreq(self): #gerer les ^ ?
	"""Return the observed frequencies of dinucleotides in the
	EvolSequence.
        
	"""
	
	c=Compte()
	c.add_seq(self,2)
	c.__idiv__(float(len(self)-1))
	return c

    def __str_freq(self): #methode pas propre
	"""Return the observed nucleotide frequencies in a simple
	string format: each value frequency on a different column.

        """
      
	l=str(self.freq()).split()
	s=""
#	s="A\tC\tG\tT\n"
	for j in string.split('A C G T'):
	    if l.count(j)==1:
		s+=l[l.index(j)+1]+'\t'
	    else:
		s+='\t0\t' 
	return(s+'\n')
   
    def __str_difreq(self): #methode pas propre
	"""Return the observed dinucleotide frequencies in a simple
	string format: each value frequency on a different column.

        """
        
	l=str(self.difreq()).split()
	s=""
#	s="AA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n"
	if len(l)%2!=0:
	    l.pop(len(l)-1)
	for j in string.split('A C G T'):
	    if l.count(j)==1: #pour virer le compte qui correspond aux bords...
		i=l.index(j)
		l.pop(i)
		l.pop(i)
	for j in string.split('AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT'):
	    if l.count(j)==1:
		s+=l[l.index(j)+1]+'\t'
	    else:
		s+='\t0\t'
	return(s+'\n')

    def __evolve_stat(self,m,d,n,nf1,nf2, **kw):
	"""Effectue des substitutions sur une séquence pendant $2
	substitutions par site et stocke les infos $3 fois, selon le
	modele de substitution encode dans la Proportion $1, puis
	calcule le freq et le difreq, et met tous ces resultats dans
	un joli fichier de nucleo $4 et de dinucleo $5"""

	nucl=""
	dinucl=""
	d=float(d)/n
	for i in range(n):
	    self.evolve(m,d,**kw)
	    nucl+=self.__str_freq()
	    dinucl+=self.__str_difreq()
	if type(nf1)==str and type(nf2)==str:
	    f=open(nf1,'w')
	    f.write(nucl)
	    f.close()
	    f=open(nf2,'w')
	    f.write(dinucl)
	    f.close()
	else:
	    print "Specifier deux noms valides de fichiers pour l'ecriture"

#####################################
################## evolution methods:

    def evolve(self,m,d, **kw): #verifier que nb subst est bien atteint
	"""Evolve the EvolSequence, according to model $1, for $2
substitutions.

Keyword argument to determine how the number of substitutions on the
sequence is approximated: [approx=string] 'Rounded' is the default.

Keyword argument to determine the algorithm used to simulate
neighbor-dependent substitutions: [algo=string] 'Berard' is the
default.

Keyword argument 'segments' is a list of all allowed
segments. Otherwise a segment is a SORTED list [beg,end], with
end position excluded. The probability of a position to be
chosen is proportional to the number of segments this position
belongs to.

WARNING: When an empty list is given, all positions are considered
allowed.

"""

        if kw.has_key("approx"):
            approx=kw["approx"]
        else:
            approx="Rounded"

        if kw.has_key("segments"):
            segments=kw["segments"]
        else:
            segments=[]

        if kw.has_key("algo"):
            algo=kw["algo"]
        else:
            algo="Berard"

        if segments==[]:
            l=len(self)
        else:
            l=0
            for s in segments:
                l+=s[1]-s[0]
        D=d*l
	anc=self.copy() # copy ancestral sequence
	n=0 #number of iterations
	s=0 #number of substitutions
	while s<D:
            i=self.pick_position(segments=segments)
	    self.substitute(m,i,algo=algo) #try to substitute the site at position i
	    n+=1
            s+=random.expovariate(m._Model__max)

	return self

    def evolve_seg(self,dmod,d, **kw): #verifier que nb subst est bien atteint
	"""Evolve the EvolSequence, according to models in dictionary
        $1, for $2 substitutions.

        The dictionary items are (deb,fin):mod where the modele mod is
        applied to the positions in range [deb:fin] (fin is excluded).
        If ranges overlap, a random model is applied on each
        substitution that occurs in the overlap. The probability of a
        position to be substituted is proportional to the number of
        model ranges this position belongs to.

        Keyword argument to determine how the number of substitutions on the
        sequence is approximated: [approx=string] 'Rounded' is the default.

        Keyword argument to determine the algorithm used to simulate
        neighbor-dependent substitutions: [algo=string] 'Berard' is the
        default.

        """

        if kw.has_key("approx"):
            approx=kw["approx"]
        else:
            approx="Rounded"

        if kw.has_key("algo"):
            algo=kw["algo"]
        else:
            algo="Berard"

        segments=[]
        l=0
        for v in dmod.keys():
            segments.append([min(v),max(v)])
            l+=max(v)-min(v)
        D=d*l
	anc=self.copy() # copy ancestral sequence
	n=0 #number of iterations
	s=0 #number of substitutions
	while s<D:
            i=self.pick_position(segments=segments) #choose a position
            lm=[]
            for deb,fin in dmod:
                if i>=deb and i<fin:
                    lm.append(dmod[deb,fin])
            if len(lm)!=0:
                m=random.choice(lm)
                self.substitute(m,i,algo=algo) #try to substitute the site at position i
                n+=1
                s+=random.expovariate(m._Model__max)

	return self

    def substitute(self,m,pos, approx='Rounded', algo='Berard'):
	"""Substitute, according to model $1, position $2 in the
	EvolSequence. This can result in no modification of the
	nucleotide at position $2.
	"""

	r=random.uniform(0,m._Model__max) # _Model__max exclu	
        if algo=='Berard':
            a=m._Model__next_all.get(self[pos],{}) # all allowed substitutions on this position
            for k in a:
                if r<k[3]:
                    # print '**************yes*************' #substitution might take place
                    # test a gauche
                    i=0
                    mg=k[0]
                    lg=len(mg)
                    md=k[1]
                    ld=len(md)
                    while i < lg:
                        if self[pos-lg+i]!=mg[i]:
                            return k[3]
                        else:
                            i+=1
                    i=0
                    while i < ld:
                        if pos+1+i >= len(self) or self[pos+1+i]!=md[i]:
                            return 0
                        else:
                            i+=1
                            
                    if self[pos]!=k[2]:
                        self[pos]=k[2]
                        return k[3]
                    else:
                        return 0
                else:
                    r-=k[3]
        else:
            raise NotImplementedError, "Only the 'Berard' algorithm has been implemented."
	return 0

    def replace(self,i,j):
	"""Replace nucleotide in position $1 by nucleotide $2.
	"""
	try:
	    self[i]=j
	except IndexError, e:
	    raise IndexError, 'Bad index.'
	return
    
#######################################################################
#######################################################################
########  miscellaneous methods

def n_subst(d,n,approx='Rounded'): # no use anymore
    """Compute the number of substitutions that will occur on a
    sequence given $1 substitutions per nucleotide and $2 the
    sequence length.
    
    Keyword argument: [approx=string] 'Rounded' is the default.
    """

    D=d*n
    if approx=='Rounded':
        try:
            eucl=divmod(D,int(D))
        except ZeroDivisionError, e:
            D=0
            return D
        if random.random() < eucl[1]: #euclidian division modulo
            D=int(D)
        else:
            D=int(D)+1
    return D
