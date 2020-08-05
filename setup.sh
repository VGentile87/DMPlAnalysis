#!/bin/bash
export DMPLA="/home/vale/DMPlAnalysis"
export LD_LIBRARY_PATH=$DMPLA/lib:${LD_LIBRARY_PATH}
export PATH=${PATH}:$DMPLA/lib
export DMPLA_INC=$DMPLA/include

alias cpstart='cp $DMPLA/start_root6.sh .'
