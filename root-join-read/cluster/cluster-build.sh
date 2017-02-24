#!/bin/bash

# 2017-02-16

# + add internalPb210
# ~ switch to join41.py


base=/home/mkauer/mc-fitting/root-join-read
cluster=$base/cluster
build=$base/build41
join=$base/join41.py

#for mcfile in data airshield internalK40 internalTh232 internalU238 internalPb210 lsveto lsvetoair pmtK40 pmtU238 pmtTh232 steel
for mcfile in data
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
