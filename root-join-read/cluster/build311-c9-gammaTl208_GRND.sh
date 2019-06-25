#!/bin/bash

source /home/mkauer/root-env.sh

export nodename=`hostname`
export pbsuser=$USER
export pbsdate=`date +\"%Y-%m-%d\"`
export workdir=/home/mkauer/mc-fitting/root-join-read/cluster

/home/mkauer/mc-fitting/root-join-read/join311.py /home/mkauer/mc-fitting/root-join-read/build311/configs/c9-gammaTl208_GRND > /home/mkauer/mc-fitting/root-join-read/cluster/log-build311-c9-gammaTl208_GRND.log 2>&1
