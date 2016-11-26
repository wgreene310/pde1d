export REL=4.0.0
export REL=4.2.0
OCT=/cygdrive/c/Octave/Octave-$REL/bin
#OCT=/cygdrive/c/Octave/Octave-4.2.0-rc2/bin
export PATH=$OCT:$PATH
#which mkoctfile
#echo $PATH
#/usr/bin/make -f Makefile_octave $*
#C:/MinGW/msys/1.0/bin/make -f Makefile_octave
$HOME/bin/make-mingw32.exe -f Makefile_octave $*
if [ -e pde1d.mex ]; then
/cygdrive/c/Progra~2/MICROS~1.0/VC/bin/amd64/dumpbin -dependents pde1d.mex
DEST=$HOME/scripts/octave/pde1d/inst/win32
/usr/bin/mv pde1d.mex $DEST
export LD_LIBRARY_PATH=$OCT
fi
