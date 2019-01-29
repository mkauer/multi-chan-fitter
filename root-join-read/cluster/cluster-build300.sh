#!/bin/bash

# 2019-01-21

#=================================================
V=300  # remember to update the build number
#=================================================

build="build$V"
join="join$V.py"
base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
configdir=$builddir/configs
joinscript=$base/$join


#for mcfile in `ls $configdir | grep SET1-V00-04-04`
#for mcfile in `ls $configdir | grep SET1-V00-04-12`
#for mcfile in `ls $configdir | grep SET2-V00-04-12`
#for mcfile in `ls $configdir | grep -v data-`
#for mcfile in c6-data-SET2-V00-04-12 c7-data-SET2-V00-04-12 c8-data-SET2-V00-04-12
#for mcfile in `ls $configdir | grep U235`
#for mcfile in `ls $configdir | grep Te121`
for mcfile in `ls $configdir | grep -E "copper|teflon|surf"`
do
    echo -e $mcfile
    file="$clustdir/$build-$mcfile.sh"
    cat > $file <<EOF
#!/bin/bash

source /home/mkauer/root-env.sh

export nodename=\`hostname\`
export pbsuser=\$USER
export pbsdate=\`date +\"%Y-%m-%d\"\`
export workdir=$clustdir

$joinscript $configdir/$mcfile > $clustdir/log-$build-$mcfile.log 2>&1
EOF
    
    #qsub -q "very_short" $file  #( 2 hours)
    #qsub -q "short" $file       #(12 hours)
    qsub -q "medium" $file      #(48 hours)
    #qsub -q "long" $file        #(96 hours)
    
done

