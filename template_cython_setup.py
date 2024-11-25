"""XZCython

Cython:
  -> interface: *.pxd           # not allows pure python `def`s
  -> wrapper, *.pyx  &  *.py    # same name `interface` will be automatically searched
After compiled:
  -> linux: *.so
  -> windows: *.pyd
Note:
    `wrapper` should be different name with `interface`,
    otherwise, signatures inside wrapper will be used instead

Compiler directives:
    https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#compiler-directives

Some important directives:
    # distutils: language = C++
    # distutils: libraries = lib1  lib2
    # distutils: include_dirs = dir1  dir2
    # distutils: sources = s.c  s.cpp
    # cython: language_level = 3        # python version, no need
    #link: https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html

Decorators:
    @cython.exceptval(check=True)
    @cython.cfunc       # create `cdef` function
    @cython.cclass      # create `cdef class`
    @cython.ccall       # create `cpdef` function
    @cython.locals      # local variables
    @cython.inline      # equivalent with C inline
    @cython.returns
    @cython.profile
    @cython.declare

Keywords:
    `cdef`      : used for internal C functions
    `cpdef`     : visible to Python

    cdef class cls
    cdef [public|inline] [int|str|bint] var
    cdef [struct|union|enum] block

    DEF key value       # work like macros in compiling time


Compile Process:
    1) generate C/C++ source codes:
        cython source.pyx   # new file: `source.c` or `source.cpp`
    
    2) compile:
        gcc -pthread -B /home/xiang/Applications/miniconda3/compiler_compat \
            -Wsign-compare -g -fwrapv -O3 -Wall -fPIC       \
            -I/home/xiang/Applications/miniconda3/include/python3.8 \
            -I/all/my/includes \
            -c source.cpp       -o source.o
        g++ -pthread -B /home/xiang/Applications/miniconda3/compiler_compat -shared \
            -L/home/xiang/Applications/miniconda3/lib       \
            -Wl,-rpath=/home/xiang/Applications/miniconda3/lib      \
            -I/all/my/includes \
            allmyobjs.o         -o final.so
"""


from setuptools import Extension, setup
from Cython.Build import cythonize

setup(
    name='obj',
    ext_modules = cythonize([Extension("queue", ["queue.pyx"])])
)


"""Files

# interface file: cqueue.pxd
# (selectively copy-and-paste of header file)
cdef extern from "queue.h":
    ctypedef struct Queue:
        pass
    Queue* queue_new()


# wrapper file: queue.pyx
# (following compiler directives should be defined)
# distutils: sources = path/to/source/queue.c
# distutils: include_dirs = path/to/source/
cimport cqueue
cdef func(int a):
    pass
cdef class Temp:
    def __cinit__(self):
        pass
    def __dealloc__(self):
        pass
"""




