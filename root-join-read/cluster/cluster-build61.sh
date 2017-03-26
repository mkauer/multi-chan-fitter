#!/bin/bash

# 2017-03-26

# ~ break down jobs by isotope
# ~ testing new variables
# ~ new primPMTid pmt cuts
# ~ switch to join/build 61
# + add internal I125
# + add internal Na22
# + add internalsurf Pb210
# ~ switch to join60.py
# + add internal Pb210
# ~ switch to join41.py

#=================================================

V=61

build="build$V"
join="join$V.py"

#=================================================


base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
joinscript=$base/$join

for mcfile in data airshield internalNa22 internalI125 internalK40 internalTh232 internalU238 internalPb210 internalsurfPb210 lsvetoK40 lsvetoU238 lsvetoTh232 lsvetoair pmtK40 pmtU238 pmtPb210 pmtTh232 steel
do
    file="$clustdir/$build-$mcfile.sh"
    cat > $file <<EOF
#!/bin/bash

source /home/mkauer/env.sh

export nodename=\`hostname\`
export pbsuser=\$USER
export pbsdate=\`date +\"%Y-%m-%d\"\`
export workdir=$clustdir

$joinscript $builddir/$mcfile > $clustdir/log-$build-$mcfile.log 2>&1
EOF
    
    #qsub -q "very_short" $file
    qsub -q "short" $file
    
done

### qsub -q
# "very_short"
# "short"
# "medium"
# "long"
