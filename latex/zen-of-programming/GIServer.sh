#!/bin/bash
if [ $1 -eq 4 ]
then 
        IP=00.00.0.000
elif [ $1 -eq 1 ]
then
        IP=00.00.0.001
else
        exit
fi

rdesktop -d domain -g 2400x1300 -u emilioberti $IP

