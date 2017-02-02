#!/bin/bash

# 2017-02-02

### qsub -q
# "very_short"
# "short"
# "medium"
# "long"

base=/home/mkauer/mc-fitting/root-join-read
cluster=$base/cluster
build=$base/build40

for mcfile in airshield  data  internalK40  internalTh232  internalU238  lsveto  lsvetoair  pmt  steel; do
    
    file="$cluster/build-$mcfile.sh"
    cat > $file <<EOF
#!/bin/bash

source /home/mkauer/env.sh

export nodename=\`hostname\`
export pbsuser=\$USER
export pbsdate=\`date +\"%Y-%m-%d\"\`
export workdir=$cluster

$base/join40.py $build/$mcfile > $cluster/log-build-$mcfile.log 2>&1
EOF
    
    #qsub -q "very_short" $file
    qsub -q "short" $file
    
done

