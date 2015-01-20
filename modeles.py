def JC(**kw):
    """ Jukes and Cantor 1969.
    Normalised such that the substitution rate on a stationnary
sequence with no-dependency is 1.

rCgT: relative rate (with C->T) of neighbour dependant
     CpG -> TpG mutations (default: 0)
rcGA: relative rate (with G->A) of neighbour dependant
     CpG -> CpA mutations (default: 0)
eta: rate (default=1)
"""

    rCgT=kw.get("rCgT",0)
    rcGA=kw.get("rcGA",0)
    eta=kw.get("eta",1)

    for s in ["rCgT","rcGA","eta"]:
        exec("if "+s+" <0: raise ValueError, 'Bad value for "+s+"'")

    eta2=eta/3.0
    
    s=""
    s+="A|C %f\n"%eta2
    s+="A|G %f\n"%eta2
    s+="A|T %f\n"%eta2
    s+="C|A %f\n"%eta2
    s+="C|G %f\n"%eta2
    s+="C|T %f\n"%eta2
    s+="G|A %f\n"%eta2
    s+="G|C %f\n"%eta2
    s+="G|T %f\n"%eta2
    s+="T|A %f\n"%eta2
    s+="T|C %f\n"%eta2
    s+="T|G %f\n"%eta2
    s+="Cg|T %f\n"%(rCgT*eta2)
    s+="cG|A %f\n"%(rCgT*eta2)

    return s

def K80(**kw):
    """ Kimura 1980.
    Normalised such that the substitution rate on a stationnary
sequence with no-dependency is 1.

kappa: transition/transversion rate.
rCgT: relative rate (with C->T) of neighbour dependant
     CpG -> TpG mutations (default: 0)
rcGA: relative rate (with G->A) of neighbour dependant
     CpG -> CpA mutations (default: 0)
eta: rate (default=1)
"""
    kappa=kw.get("kappa",1)
    rCgT=kw.get("rCgT",0)
    rcGA=kw.get("rcGA",0)
    eta=kw.get("eta",1)

    for s in ["kappa","rCgT","rcGA","eta"]:
        exec("if "+s+" <0: raise ValueError, 'Bad value for "+s+"'")
    
    eta2=eta/(kappa+2.0)

    s=""
    s+="A|C %f\n"%eta2
    s+="A|G %f\n"%(kappa*eta2)
    s+="A|T %f\n"%eta2
    s+="C|A %f\n"%eta2
    s+="C|G %f\n"%eta2
    s+="C|T %f\n"%(kappa*eta2)
    s+="G|A %f\n"%(kappa*eta2)
    s+="G|C %f\n"%eta2
    s+="G|T %f\n"%eta2
    s+="T|A %f\n"%eta2
    s+="T|C %f\n"%(kappa*eta2)
    s+="T|G %f\n"%eta2
    if rCgT!=0:
        s+="Cg|T %f\n"%(rCgT*kappa*eta2)
    if rcGA!=0:
        s+="cG|A %f\n"%(rcGA*kappa*eta2)
    
    return s


def T92(**kw):
    """ Tamura 1992.
    Normalised such that the substitution rate on a stationnary
sequence with no-dependency is 1.

kappa: transition/transversion rate (default: 1).
theta: equilibrium frequency of GC (default: 0.5).
rCgT: relative rate (with C->T) of neighbour dependant
     CpG -> TpG mutations (default: 0)
rcGA: relative rate (with G->A) of neighbour dependant
     CpG -> CpA mutations (default: 0)
eta: relative rate with no-dependant model(default=1)
"""
    rCgT=kw.get("rCgT",0)
    rcGA=kw.get("rcGA",0)
    kappa=kw.get("kappa",1)
    theta=kw.get("theta",0.5)
    eta=kw.get("eta",1)
    
    for s in ["kappa","theta", "rCgT","rcGA", "eta"]:
        exec("if "+s+" <0: raise ValueError, 'Bad value for "+s+"'")
        
    for s in ["theta"]:
        exec("if "+s+" >1: raise ValueError, 'Bad value for "+s+"'")

    eta2=eta/(1.0+2*theta*kappa-2*theta*theta*kappa)

    s=""
    s+="A|C %f\n"%(theta*eta2)
    s+="A|G %f\n"%(theta*kappa*eta2)
    s+="A|T %f\n"%((1-theta)*eta2)
    s+="C|A %f\n"%((1-theta)*eta2)
    s+="C|G %f\n"%(theta*eta2)
    s+="C|T %f\n"%((1-theta)*kappa*eta2)
    s+="G|A %f\n"%((1-theta)*kappa*eta2)
    s+="G|C %f\n"%(theta*eta2)
    s+="G|T %f\n"%((1-theta)*eta2)
    s+="T|A %f\n"%((1-theta)*eta2)
    s+="T|C %f\n"%(theta*kappa*eta2)
    s+="T|G %f\n"%(theta*eta2)
    if rCgT!=0:
        s+="Cg|T %f\n"%(rCgT*(1-theta)*kappa*eta2)
    if rcGA!=0:
        s+="cG|A %f\n"%(rcGA*(1-theta)*kappa*eta2)

    return s


def HKY85(**kw):

    """ Hasegawa Kishino Yano 1985
    Normalised such that the substitution rate on a stationnary
sequence with no-dependency is 1.

kappa: transition/transversion rate
theta: equilibrium frequency of GC.
theta1: ratio of equilibrium frequency of A on equilibrium
    frequency of AT
theta2: ratio of equilibrium frequency of G on equilibrium
    frequency of CG

rCgT: relative rate (with C->T) of neighbour dependant
    CpG -> TpG mutations (default: 0)
rcGA: relative rate (with G->A) of neighbour dependant
    CpG -> CpA mutations (default: 0)
eta: rate (default=1)
"""
    kappa=kw.get("kappa",1)
    theta=kw.get("theta",0.5)
    theta1=kw.get("theta1",0.5)
    theta2=kw.get("theta2",0.5)
    rCgT=kw.get("rCgT",0)
    rcGA=kw.get("rcGA",0)
    eta=kw.get("eta",1)
    
    for s in ["kappa","theta","theta1","theta2","rCgT","rcGA","eta"]:
        exec("if "+s+" <0: raise ValueError, 'Bad value for "+s+"'")

    for s in ["theta","theta1","theta2"]:
        exec("if "+s+" >1: raise ValueError, 'Bad value for "+s+"'")

    piA=theta1*(1-theta)
    piC=(1-theta2)*theta
    piG=theta2*theta
    piT=(1-theta1)*(1-theta)

    eta2=eta/(2.0*(piA*piC+piC*piG+piA*piT+piG*piT+kappa*(piC*piT+piA*piG)))

    s=""
    s+="A|C %f\n"%(piC*eta2)
    s+="A|G %f\n"%(kappa*piG*eta2)
    s+="A|T %f\n"%(piT*eta2)
    s+="C|A %f\n"%(piA*eta2)
    s+="C|G %f\n"%(piG*eta2)
    s+="C|T %f\n"%(kappa*piT*eta2)
    s+="G|A %f\n"%(kappa*piA*eta2)
    s+="G|C %f\n"%(piC*eta2)
    s+="G|T %f\n"%(piT*eta2)
    s+="T|A %f\n"%(piA*eta2)
    s+="T|C %f\n"%(kappa*piC*eta2)
    s+="T|G %f\n"%(piG*eta2)
    if rCgT!=0:    
        s+="Cg|T %f\n"%(rCgT*kappa*piT*eta2)
    if rcGA!=0:
        s+="cG|A %f\n"%(rcGA*kappa*piT*eta2)

    return s



def TN93(**kw):
    """ Tamura and Nei 1993.
    Normalised such that the substitution rate on a stationnary
sequence with no-dependency is 1.

kappa1: transition/transversion for A<->G.
kappa2: transition/transversion for C<->T.
theta: equilibrium frequency of GC.
theta1: ratio of equilibrium frequency of A on equilibrium
     frequency of AT
theta2: ratio of equilibrium frequency of G on equilibrium
     frequency of CG
    
rCT: relative rate (with C->T) of neighbour dependant
     CpG -> TpG mutations.
eta: rate (default=1)
"""

    rCgT=kw.get("rCgT",0)
    rcGA=kw.get("rcGA",0)
    kappa1=kw.get("kappa1",1)
    kappa2=kw.get("kappa2",1)
    theta=kw.get("theta",0.5)
    theta1=kw.get("theta1",0.5)
    theta2=kw.get("theta2",0.5)
    eta=kw.get("eta",1)

    for s in ["kappa1","kappa2","theta","theta1","theta2","rCgT","rcGA","eta"]:
        exec("if "+s+" <0: raise ValueError, 'Bad value for "+s+"'")

    for s in ["theta","theta1","theta2"]:
        exec("if "+s+" >1: raise ValueError, 'Bad value for "+s+"'")

    piA=theta1*(1-theta)
    piC=(1-theta2)*theta
    piG=theta2*theta
    piT=(1-theta1)*(1-theta)

    eta2=eta/(2.0*(piA*piC+piC*piG+piA*piT+piG*piT+kappa2*piC*piT+kappa1*piA*piG))

    s=""
    s+="A|C %f\n"%(piC*eta2)
    s+="A|G %f\n"%(kappa1*piG*eta2)
    s+="A|T %f\n"%(piT*eta2)
    s+="C|A %f\n"%(piA*eta2)
    s+="C|G %f\n"%(piG*eta2)
    s+="C|T %f\n"%(kappa2*piT*eta2)
    s+="G|A %f\n"%(kappa1*piA*eta2)
    s+="G|C %f\n"%(piC*eta2)
    s+="G|T %f\n"%(piT*eta2)
    s+="T|A %f\n"%(piA*eta2)
    s+="T|C %f\n"%(kappa2*piC*eta2)
    s+="T|G %f\n"%(piG*eta2)
    if rCgT!=0:
        s+="Cg|T %f\n"%(rCgT*kappa2*piT*eta2)
    if rcGA!=0:
        s+="cG|A %f\n"%(rcGA*kappa2*piT*eta2)

    return s



def F84(*kw):
    
    """ Felsenstein 1984.
    Normalised such that the substitution rate on a stationnary
sequence with no-dependency is 1.

kappa:  ratio transition/transversion
theta:  equilibrium frequency of GC.
theta1: ratio of equilibrium frequency of A on equilibrium
     frequency of AT
theta2: ratio of equilibrium frequency of G on equilibrium
     frequency of CG
    
rCT: relative rate (with C->T) of neighbour dependant
     CpG -> TpG mutations.
eta: rate (default=1)
"""
    kappa=kw.get("kappa",1)
    theta=kw.get("theta",0.5)
    theta1=kw.get("theta1",0.5)
    theta2=kw.get("theta2",0.5)
    rCgT=kw.get("rCgT",0)
    rcGA=kw.get("rcGA",0)
    eta=kw.get("eta",1)

    for s in ["kappa","theta","theta1","theta2","rCgT","rcGA","eta"]:
        exec("if "+s+" <0: raise ValueError, 'Bad value for "+s+"'")

    for s in ["theta","theta1","theta2"]:
        exec("if "+s+" >1: raise ValueError, 'Bad value for "+s+"'")

    piA=theta1*(1-theta)
    piC=(1-theta2)*theta
    piG=theta2*theta
    piT=(1-theta1)*(1-theta)
    piR=piA+piG
    piY=piC+piT
    
    eta2=eta/(2.0*kappa*(piC*piT/(piC+piT)+piA*piG/(piA+piG))-piC*piC-piG*piG-piT*piT-piA*piA+1)

    s=""
    s+="A|C %f\n"%(piC*eta2)
    s+="A|G %f\n"%((1+kappa/piR)*piG*eta2)
    s+="A|T %f\n"%(piT*eta2)
    s+="C|A %f\n"%(piA*eta2)
    s+="C|G %f\n"%(piG*eta2)
    s+="C|T %f\n"%((1+kappa/piY)*piT*eta2)
    s+="G|A %f\n"%((1+kappa/piR)*piA*eta2)
    s+="G|C %f\n"%(piC*eta2)
    s+="G|T %f\n"%(piT*eta2)
    s+="T|A %f\n"%(piA*eta2)
    s+="T|C %f\n"%((1+kappa/piY)*piC*eta2)
    s+="T|G %f\n"%(piG*eta2)
    if rCgT!=0:
        s+="Cg|T %f\n"%(rCgT*(1+kappa/piY)*piT*eta2)
    if rcGA!=0:
        s+="cG|A %f\n"%(rcGA*(1+kappa/piY)*piT*eta2)

    return s



def GTR(**kw):
    """ General Time-Reversible substitution model,
parametrization of Yang 1994.
    Normalised such that the substitution rate on a stationnary
sequence with no-dependency is 1.

In proportion to the A<->G rate, the arguments are (default 1):

a: C<->T
b: A<->T
c: G<->T
d: A<->C
e: C<->G

theta: equilibrium frequency of GC.
theta1: ratio of equilibrium frequency of A on equilibrium
     frequency of AT
theta2: ratio of equilibrium frequency of G on equilibrium
     frequency of CG
    
rCgT: relative rate (with C->T) of neighbour dependant
     CpG -> TpG mutations (default: 0)
rcGA: relative rate (with G->A) of neighbour dependant
     CpG -> CpA mutations (default: 0)
eta: rate (default=1)
"""

    a=kw.get("a",1)
    b=kw.get("b",1)
    c=kw.get("c",1)
    d=kw.get("d",1)
    e=kw.get("e",1)
    theta=kw.get("theta",0.5)
    theta1=kw.get("theta1",0.5)
    theta2=kw.get("theta2",0.5)
    rCgT=kw.get("rCgT",0)
    rcGA=kw.get("rcGA",0)
    eta=kw.get("eta",1)

    for s in ["a", "b", "c", "d", "e", "theta","theta1","theta2","rCgT", "rcGA", "eta"]:
        exec("if "+s+" <0: raise ValueError, 'Bad value for "+s+"'")

    for s in ["theta","theta1","theta2"]:
        exec("if "+s+" >1: raise ValueError, 'Bad value for "+s+"'")

    piA=theta1*(1-theta)
    piC=(1-theta2)*theta
    piG=theta2*theta
    piT=(1-theta1)*(1-theta)
    piR=piA+piG
    piY=piC+piT

    eta2=eta/(2.0*(a*piC*piT+b*piA*piT+c*piC*piG+d*piA*piC+e*piC*piG+piA*piG))
    
    s=""
    s+="A|C %f\n"%(d*piC*eta2)
    s+="A|G %f\n"%(piG*eta2)
    s+="A|T %f\n"%(b*piT*eta2)
    s+="C|A %f\n"%(d*piA*eta2)
    s+="C|G %f\n"%(e*piG*eta2)
    s+="C|T %f\n"%(a*piT*eta2)
    s+="G|A %f\n"%(piA*eta2)
    s+="G|C %f\n"%(e*piC*eta2)
    s+="G|T %f\n"%(c*piT*eta2)
    s+="T|A %f\n"%(b*piA*eta2)
    s+="T|C %f\n"%(a*piC*eta2)
    s+="T|G %f\n"%(c*piG*eta2)
    if rCgT!=0:
        s+="Cg|T %f\n"%(rCgT*a*piT*eta2)
    if rcGA!=0:
        s+="cG|A %f\n"%(rcGA*a*piT*eta2)

    return s

