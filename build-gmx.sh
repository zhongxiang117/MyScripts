
# be aware: those configs have to be run before installation starts!!!
unset CC  CC_FOR_BUILD  CFLAGS
unset CMAKE_ARGS  CMAKE_PREFIX_PATH
unset CPP  CPPFLAGS  CXX  CXXFILT  CXXFLAGS  CXX_FOR_BUILD
unset GCC  GCC_AR  GCC_NM  GCC_RANLIB  GXX
unset LD  LDFLAGS  LD_GOLD
unset NVCC_PREPEND_FLAGS
unset DESTDIR  CPATH

# note: version CUDA 12 -> GCC 12
# GMX 2025 require: CUDA >=12, thus, for CUDA=11, GMX 2024.5 is the last version to use


echo $NVCC_PREPEND_FLAGS
export NVCC_PREPEND_FLAGS=-I/home/xzhong/Applications/miniconda3/envs/utils3/targets/x86_64-linux/include

if [[ -f CMakeCache.txt ]]; then
    echo "exist: CMakeCache.txt, delete? y/n"
    read tmp
    if [[ $tmp != y ]]; then exit; fi
    #rm -f CMakeCache.txt
fi

# for Plumed and MPI
# Plumed: use v2.9.0 but not v2.9.1, reason: https://github.com/plumed/plumed2/issues/981
export LDFLAGS="-L/home/xzhong/Applicatoins/build-plumed-v2.9/lib -L/home/xzhong/Applications/miniconda3/envs/utils2/lib"
export CMAKE_INCLUDE_PATH=/home/xzhong/Applications/build-plumed-v2.9/include

# for 'cuda_runtime.h' not found
export CPATH=/home/xzhong/Applications/miniconda3/envs/utils3/targets/x86_64-linux/include

cmake ../ \
    -DGMX_GPU=CUDA \
    -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON \
    -DCUDA_TOOLKIT_INCLUDE=/home/xzhong/Applications/miniconda3/envs/utils3/targets/x86_64-linux/include/ \
    --debug-output \


make -j 20
export DESTDIR=where-to-install; make install

# when make failed, to show full commands by examining
#make VERBOSE=1
#make -n


# When error:
#make[2]: *** No rule to make target `/usr/local/lib/libplumed.so', needed by `bin/gmx'.  Stop.
#make[1]: *** [src/programs/CMakeFiles/gmx.dir/all] Error 2
#make: *** [all] Error 2

#> from `make[1]`, manually set target file in "src/programs/CMakeFiles/gmx.dir/build.make"


# /home/xiang/miniconda3/envs/xutils/x86_64-conda-linux-gnu/bin/ld: cannot find -lplumed: No such file or directory
#collect2: error: ld returned 1 exit status
#make[2]: *** [bin/gmx] Error 1
#make[1]: *** [src/programs/CMakeFiles/gmx.dir/all] Error 2
#make: *** [all] Error 

#> similarly, manually set "-L" in "src/programs/CMakeFiles/gmx.dir/link.txt"


