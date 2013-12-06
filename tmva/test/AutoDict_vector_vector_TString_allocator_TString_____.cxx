#include "vector"
#include "vector"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<vector<TString,allocator<TString> > >+;
#pragma link C++ class vector<vector<TString,allocator<TString> > >::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<vector<TString,allocator<TString> > >::iterator;
#pragma link C++ operators vector<vector<TString,allocator<TString> > >::const_iterator;
#pragma link C++ operators vector<vector<TString,allocator<TString> > >::reverse_iterator;
#endif
#endif
