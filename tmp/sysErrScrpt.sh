#!/bin/sh
#nJob1=$1
######################################################
#
# This script is for the systematic analysis of flow analysis for Run 18 7.2 GeV
# Flow analysis
# author: Ding Chen
# date: 9/15/2020
#
######################################################
nJob1=1
maxJob1=$1
# how many systematic cuts you have: 13 for flow analysis in this script

nJob2=1
maxJob2=$2
# how many variations for each systematic cut: 1 or 2

nJob3=$3
nNumber=0 # total # of submits

# which iteration of the flow calculation: 1-3

# submit scheduler jobs in batch
DEBUG=true
while [ $nJob1 -le $maxJob1 ]
do
    echo $nJob1
    nJob2=1
    while [ $nJob2 -le $maxJob2 ]
    do
        echo $nJob2
        if ${DEBUG}; then
        # Print the command
        echo "star-submit-template -template submitPicoDstJobs_v2_sysErr.xml -entities cutID=$nJob1,varID=$nJob2,iterID=$nJob3"
        # Execute the command
        star-submit-template -template submitPicoDstJobs_v2_sysErr.xml -entities cutID=$nJob1,varID=$nJob2,iterID=$nJob3
        ((nNumber++))
        echo "total submits: $nNumber"
        fi
        if(($nJob1 == 12 && $nJob2 ==1)); then
        break # cutID = 13: Dip Angle - 2 variations
        fi
        ((nJob2++))
    done
    ((nJob1++))
done
#echo $nJob1
