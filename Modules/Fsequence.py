# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Fsequence

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class Fsequence(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Fsequence, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Fsequence, name)
    def __repr__(self):
        return "<C Fsequence instance at %s>" % (self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Fsequence, 'this', _Fsequence.new_Fsequence(*args))
        _swig_setattr(self, Fsequence, 'thisown', 1)
    def __del__(self, destroy=_Fsequence.delete_Fsequence):
        try:
            if self.thisown: destroy(self)
        except: pass
    def melange(*args): return _Fsequence.Fsequence_melange(*args)
    def termine(*args): return _Fsequence.Fsequence_termine(*args)
    def vtaille(*args): return _Fsequence.Fsequence_vtaille(*args)
    def __getitem__(*args): return _Fsequence.Fsequence___getitem__(*args)
    def __setitem__(*args): return _Fsequence.Fsequence___setitem__(*args)
    def __str__(*args): return _Fsequence.Fsequence___str__(*args)
    def __len__(*args): return _Fsequence.Fsequence___len__(*args)
    def recup_rel(*args): return _Fsequence.Fsequence_recup_rel(*args)

class FsequencePtr(Fsequence):
    def __init__(self, this):
        _swig_setattr(self, Fsequence, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Fsequence, 'thisown', 0)
        _swig_setattr(self, Fsequence,self.__class__,Fsequence)
_Fsequence.Fsequence_swigregister(FsequencePtr)


