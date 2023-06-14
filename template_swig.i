%define DOCSTRING
"##############################################################################\n\
# Example to show how to use SWIG                                              \n\
##############################################################################\n"
%enddef

// use `swig -swiglib` to show standard library path

// two files will be created: 1) MOD.py (interface), 2) MOD_wrap.cxx
// important: its sona file should be named as _MOD.so
%module(docstring=DOCSTRING, package="test") MOD

// all things including comments will be copied as-is at the beginning of MOD_wrap.cxx
%begin %{
#define SWIG_PYTHON_2_UNICODE
//#define SWIG_FILE_WITH_INIT
%}


// Standard blocks help on most C/C++ source codes
//#############################################################################

// Add standard C++ library
%include "std_array.i"
%include "std_list.i"
%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"
%include "std_pair.i"
%include "typemaps.i"

%rename(inc) *::operator++;
%ignore *::operator=;
%ignore *::operator[];

// Help SWIG to understand some special types, like list of strings
namespace std {
    %template(IntVector) vector<int>;
    %template(DoubleVector) vector<double>;
    %template(DoubleVectorVector) vector<vector<double>>;
    %template(StringVector) vector<string>;
    %template(ConstCharVector) vector<const char*>;
}

typedef int size_t;
typedef int sz;
typedef double fl;

//#############################################################################


// ignore the methods of vector operation in target language
// note: ## in C/C++ works as the concatenation
//%define VECTORTEMPLATE_WRAP(vectorname, T)
//%feature("ignore") vector<T>::append;
//%feature("ignore") vector<T>::assign;
//%feature("ignore") vector<T>::back;
//%feature("ignore") vector<T>::begin;
//%feature("ignore") vector<T>::capacity;
//%feature("ignore") vector<T>::empty;
//%feature("ignore") vector<T>::end;
//%feature("ignore") vector<T>::erase;
//%feature("ignore") vector<T>::front;
//%feature("ignore") vector<T>::get_allocator;
//%feature("ignore") vector<T>::insert;
//%feature("ignore") vector<T>::pop;
//%feature("ignore") vector<T>::pop_back;
//%feature("ignore") vector<T>::push_back;
//%feature("ignore") vector<T>::rbegin;
//%feature("ignore") vector<T>::rend;
//%feature("ignore") vector<T>::reserve;
//%feature("ignore") vector<T>::resize;
//%feature("ignore") vector<T>::size;
//%feature("ignore") vector<T>::swap;
//%template(vector ## vectorname) vector<T>;
//%enddef


// codes will be inserted at the beginning of the python wrapper
%pythonbegin %{
import sys
if sys.platform.find("linux") != -1:
    dlflags = sys.getdlopenflags()
    import ctypes
    sys.setdlopenflags(dlflags | ctypes.RTLD_GLOBAL)
%}

// codes will be inserted after the SWIG initialization
%pythoncode %{
if sys.platform.find("linux") != -1:
    sys.setdlopenflags(dlflags)
%}


// entries inside `%{` `%}` will be copied as-is
// only used for declarations, will not be parsed
%{
#include "headers.h"
%}

%include "exception.i"
%exception {
    try {
        $action
    } catch (const error& e) {
        SWIG_exception(SWIG_TypeError, e.what());
    }
}


// rename long/ambiguous oldFunc to new simple consice newFunc
//%rename(newFunc) oldFunc;


%{
#include "specific.h"
%}


// %ignore someFunc;        // only funcname is enough


// include only once, parsing declarations entry in here
%include "specific.h"


/*
## Steps to SWIG:

## new file: interface.py, interface_wrap.cxx
swig -python -c++ -small -naturalvar -fastdispatch -shadow -py3 \
    -Isrc/lib \
    interface.i

test $? == 0 || { echo "Fatal: in swig"; exit 1; }


## src dir should be searched firstly to avoid name collisions
g++ -pthread -shared -g -O3 -fPIC -std=c++11 \
    -Isrc/lib   \
    -I/home/xiang/Applications/miniconda3/include/python3.8     \
    -I/home/xiang/Applications/miniconda3/include  \
    -c interface_wrap.cxx -o interface_wrap.o

test $? == 0 || { echo "Fatal: in object"; exit 1; }


## option: rpath: add a directory to the runtime library search path
g++ -pthread -shared -Wl,-rpath=/home/xiang/Applications/miniconda3/lib \
    -L/home/xiang/Applications/miniconda3/lib   \
    src/lib/*.o    \
    interface_wrap.o    \
    -o _interface_wrap.so

test $? == 0 || { echo "Fatal: in link"; exit 1; }

*/




