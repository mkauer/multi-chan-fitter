#!/bin/bash

# 2017-05-23

# ~ BIG bug fixed - wasn't selecting volumeCut for non-internal mc!
# + added support for including other pmts ie 'extpmt'
# + split internal u238 and th232
# + split pmt u238 and th232
# ~ switch to join62.py
# + add internalTe121m, internalTe123m, internalTe125m, internalTe127m
# + add data1544
# ~ mv data to data1546
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

V=63

build="build$V"
join="join$V.py"

#=================================================


base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
joinscript=$base/$join

for mcfile in extpmtU238 extpmtTh232 extpmtK40
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
