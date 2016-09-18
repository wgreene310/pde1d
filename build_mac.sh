export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
#./run-octave
octdir=/Applications/Octave.app/Contents/Resources/usr/bin
export PATH=$octdir/:$PATH
#which octave
make -f Makefile_octave $*
if [ -e pde1d.mex ]; then
otool -L pde1d.mex 
mv pde1d.mex $HOME/win7/scripts/octave/pde1d/inst/macos32
fi
