#!/bin/bash

# 2017-03-24

# ~ new primPMTid pmt cuts
# ~ switch to join/build 61
# + add internal I125
# + add internal Na22
# + add internalsurf Pb210
# ~ switch to join60.py
# + add internal Pb210
# ~ switch to join41.py


base=/home/mkauer/mc-fitting/root-join-read
cluster=$base/cluster
build=$base/build61
join=$base/join61.py

#for mcfile in data airshield internalNa22 internalI125 internalK40 internalTh232 internalU238 internalPb210 internalsurfPb210 lsveto lsvetoair pmtK40 pmtU238 pmtTh232 steel
for mcfile in pmtU238 pmtTh232
do
    file="$cluster/build-$mcfile.sh"
    cat > $file <<EOF
#!/bin/bash

source /home/mkauer/env.sh

export nodename=\`hostname\`
export pbsuser=\$USER
export pbsdate=\`date +\"%Y-%m-%d\"\`
export workdir=$cluster

$join $build/$mcfile > $cluster/log-build-$mcfile.log 2>&1
EOF
    
    #qsub -q "very_short" $file
    qsub -q "short" $file
    
done

### qsub -q
# "very_short"
# "short"
# "medium"
# "long"
