#!/bin/sh
# usage: source setup_trees_80X.csh <rel53X> <rel80X>

if [[ "X$1" == "X" ]]; then 
    echo "usage: source setup_trees_80X.sh <rel53X> <rel80X>"
    return
fi

if [[ "X$2" == "X" ]]; then 
    echo "usage: source setup_trees_80X.sh <rel53X> <rel80X>"
    return
fi



export rel53X=$1
export rel80X=$2

cd $rel53X 
eval `scramv1 runtime -sh`
export pythonPath53X=$CMSSW_BASE/python/
echo "Now in: " $rel53X

cd $rel80X 
eval `scramv1 runtime -sh`
export PYTHONPATH=$pythonPath53X\:$PYTHONPATH
echo "Now in: " $rel80X

cd $rel53X
echo "back to " $rel53X
