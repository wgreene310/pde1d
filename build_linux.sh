export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
#./run-octave
export PATH=/opt/octave-4.0.0/bin/:$PATH
#which octave
make -f Makefile_octave $*
DEST=$HOME/win7/scripts/octave/pde1d/inst/linux32
if [ -e pde1d.mex ]; then
mv pde1d.mex $DEST
ldd $DEST/pde1d.mex
fi
