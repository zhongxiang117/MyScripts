"""XZSetuptools"""

from Cython.Build import cythonize

from setuptools import Extension, setup
from distutils.command.build import build
from setuptools.command.build_ext import build_ext
import glob
import sys
import os


class xbuild(build):
    def ensure_finalized(self):
        # build_platlib: 'build/lib.linux-x86_64-3.8'
        # build_purelib: 'build/lib.linux-x86_64-3.8'
        # build_temp: 'build/temp.linux-x86_64-3.8'
        # build_scripts: 'build/scripts-3.8'
        super().ensure_finalized()


class xbuild_ext(build_ext):
    """Execution Sequence:
    
    build_ext.run
        self.compiler
        distutils.sysconfig.customize_compiler      # not make sense
        self.build_extensions()             # parallel or serial
            self.build_extension(ext)
                self.swig_sources(sources,extension)
                    -> self.compiler.spawn(
                objects = self.compiler.compile(
                    -> self.compiler._setup_compile(
                    -> self.compiler._get_cc_args(
                    -> self.compiler._compile(      # hack will happen in here
                self.compiler.link_shared_object(objects, ...
                    -> self.link                    # differentiate Linux or Windows compiler
                    -> setuptools._distutils.unixccompiler.UnixCCompiler
                    -> self.compiler.spawn(
                        -> distutils.spawn.spawn(
    """
    def build_extension(self, ext):
        # ext: setuptools.exetension.Extension('user_defined_obj')
        #
        #self.build_lib
        #self.build_temp
        #self.include_dirs
        #self.library_dirs
        #self.compiler                              # used by every extension
        self.orig_compile = self.compiler.compile
        self.orig_compile_ = self.compiler._compile
        self.compiler.compile = self.xzcompile      # hack
        self.compiler._compile = self.xz_compile    # hack
        #sources = ext.sources
        super().build_extension(ext)                # call `self.compiler.compile`

    def xzcompile(
            self, sources, output_dir=None, macros=None, include_dirs=None,
            debug=0, extra_preargs=None, extra_postargs=None, depends=None
        ):
        # `self._compile` will be called for every `sources`
        #
        # path of final output objects will be `os.path.join(self.build_temp,output_dir,${obj})`
        objects = self.orig_compile(
            sources, output_dir, macros, include_dirs, debug,
            extra_preargs, extra_postargs, depends
        )
        # `objects`` will be finally compiled and linked
        return objects

    def xz_compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
        # obj: ${build_temp} ${output_dir} ${obj}       # path to `-o`
        # src: path/to/source/code                      # path to `-c`
        # ext: '.cc' or '.c' or '.i'        # source code extension
        # cc_args: ['-DDEBUG=1', ]          # list of macros and includes
        # extra_postargs: []                # list of arguments
        # pp_opts: []                       # preinfo for `cc_args`, not used
        #
        # if `object` exists and its `mtime` newer than source, skip
        if os.path.isfile(obj):
            mtime = os.path.getmtime(obj)
            ftime = os.path.getmtime(src)
            if mtime > ftime:
                return
        self.orig_compile_(obj, src, ext, cc_args, extra_postargs, pp_opts)


"""
Extension:
    include_dirs: list of includes dirs, in source: should be written as `#include <myincludes>`
    define_macros, undef_macros: list of pair of macros, in form: `(name, value)`
        -> they work as: `#define name value` and `#undef name value`
    libraries: list of stardard libraries to search
    library_dirs: list of library
    runtime_library_dirs: for shared dynamically loaded libraries
    extra_compile_args, extra_link_args: arguments for compiler and linker
"""
ext = Extension(
    name='xz.ext',      # python dotted name, will finally be `xz/ext.so`
    sources=['../path/relative/to/setup.py/source/file'],   # C, C++, or SWIG(i) files
    include_dirs=[],
    define_macros=[],
    undef_macros=[],
    library_dirs=[],
    libraries=[],       # stardard libraries, e.g., `c` will be interpreted as `libc.so`
    runtime_library_dirs=[],
    extra_objects=[],
    extra_compile_args=[],
    extra_link_args=[],
    depends=[],
    optional=False,     # `False` means if build fails, exit: the extension is not optional
    swig_opts=[],       # works when `SWIG` source is detected
)


setup(
    name='extpkg',
    # which contains a list of `Extension`s, called by command `build_ext`
    ext_modules = cythonize([ext,]),
    # all entries can be found in `setuptools.command.~`
    cmdclass={'build':xbuild, 'build_ext':xbuild_ext,},
)



