#!/bin/bash

# 2018-06-28

# ~ use the medium queue
# ~ regenerate everything!!!

#=================================================

V=100

build="build$V"
join="join$V.py"

#=================================================


base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
configdir=$builddir/configs
joinscript=$base/$join

#for mcfile in `ls $configdir`
#for mcfile in c9-dataset1
#for mcfile in `ls $configdir | grep -i pb210_grnd`
for mcfile in `ls $configdir | grep -v dataset`
do
    echo -e $mcfile
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
    #qsub -q "short" $file
    qsub -q "medium" $file
    
done


### qsub -q
# "very_short"
# "short"  - 12 hours
# "medium" - 48 hours
# "long"

