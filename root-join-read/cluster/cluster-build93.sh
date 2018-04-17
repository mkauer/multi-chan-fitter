#!/bin/bash

# 2018-04-17

# + add newer copper Co60
# + add new Cd105 and Co60 backgrounds
# ~ check different data offsets in c6,c7
# ~ new v00-04-04 dataset
# ~ rebuild C3 PMT sim with other pmts added in
# ~ modified for just reprocessing of C8 with new resolution
# + new data set V00-04-03
# + add steel Th232 and U238
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

V=93

build="build$V"
join="join$V.py"

#=================================================


base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
configdir=$builddir/configs
joinscript=$base/$join

#for mcfile in `ls $configdir | grep -v -E "root|dataset" | grep internal`
#for mcfile in `ls $configdir | grep dataset1`
#for mcfile in `ls $configdir | grep internal`
#for mcfile in innersteelU238_GRND innersteelTh232_GRND
#for mcfile in c6-dataset1 c7-dataset1
#for mcfile in naisurf10umPb210_GRND teflonsurf2umPb210_GRND teflonbulkPb210_GRND
for mcfile in copperCo60_GRND
#for mcfile in `ls $configdir`
do
    file="$clustdir/$build-$mcfile.sh"
    cat > $file <<EOF
#!/bin/bash

source /home/mkauer/env.sh

export nodename=\`hostname\`
export pbsuser=\$USER
export pbsdate=\`date +\"%Y-%m-%d\"\`
export workdir=$clustdir

$joinscript $builddir/configs/$mcfile > $clustdir/log-$build-$mcfile.log 2>&1
EOF
    
    #qsub -q "very_short" $file
    qsub -q "short" $file
    #qsub -q "medium" $file
    
done


### qsub -q
# "very_short"
# "short" - 12 hours
# "medium"
# "long"
