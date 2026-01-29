
# brew instal gfortran fftw opencoarrays

# export ARCH_FLAGS="-arch arm64"
# export FC="gfortran $ARCH_FLAGS"
export FC="gfortran"
export XFLAG_NO_OMP='-O0  -g -fcheck=all -fbacktrace -fno-common -O3 -cpp -fcoarray=single' #  -fno-common'
# export XFLAG='-O0  -g -fcheck=all -fbacktrace -fno-common -O3 -cpp -fopenmp ' #-fcoarray=single' #  -fno-common '
# export XFLAG='-O0  -g  -fno-common -O3 -cpp -fopenmp '
export XFLAG='-O0  -g  -fno-common -O3 -cpp -fopenmp  -fbounds-check'

# export XFLAG_NO_OMP='-O3 -cpp -fopenmp -fcoarray=single' # -mcmodel=large'
# export XFLAG='-O3 -cpp -fopenmp -fcoarray=single' # -mcmodel=large'
export OFLAG_NO_OMP=$XFLAG_NO_OMP' -c'
export OFLAG=$XFLAG' -c'
# export FFTFLAG='-I/opt/homebrew/opt/fftw/include/ -lfftw3f_omp -lfftw3f -L/opt/homebrew/opt/fftw/lib'
export FFTFLAG='-I/opt/homebrew/opt/fftw/include/ -lfftw3f_omp -lfftw3f -L/opt/homebrew/opt/fftw/lib  -Wl,-rpath,/opt/homebrew/opt/fftw/lib'

export OMP_STACKSIZE=32000M
#export KMP_STACKSIZE=16000M
export OMP_NUM_THREADS=32
#export OMP_THREAD_LIMIT=4
export FOR_COARRAY_NUM_IMAGES=1
ulimit
# ulimit -s unlimited
ulimit -c unlimited  # 允许写 core dump
ulimit -v unlimited  # 内存无限制