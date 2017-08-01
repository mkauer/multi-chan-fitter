#!/bin/bash

# 2017-07-12

# + add Cd109 and H3
# ~ fixed issue with Te names not showing the 'm'
# + add LIST of files to cat and process
# ~ split data into crystal numbers
# ~ changed default histogram parameters
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

V=70

build="build$V"
join="join$V.py"

#=================================================


base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
joinscript=$base/$join

#for mcfile in `ls $builddir`
#for mcfile in dataset1
#for num in 1 2 3 4 5 6 7 8
#for mcfile in `cat TeXXXm`
#for mcfile in internalCd109_GRND internalH3_GRND
#for mcfile in pmtU238_Ra226 pmtRa226_Pb210 extpmtU238_Ra226 extpmtRa226_Pb210
#for mcfile in `cat LSVETO`
#for mcfile in pmtPb210_GRND lsvetoPb210_GRND
for mcfile in c7-dataset1
do
    #mcfile="c$num-dataset1"
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
