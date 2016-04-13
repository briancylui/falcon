export DELPHES=$PWD/delphes
export PATH=$PWD/bin:$PATH
export FASTJET=$PWD/delphes/external
export DYLD_LIBRARY_PATH=$FASTJET/lib
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH:$DELPHES
export PYTHONPATH=$PWD/python:$PYTHONPATH:$DELPHES/python



