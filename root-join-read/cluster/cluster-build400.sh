#!/bin/bash

# 2019-09-03

#=================================================
V=400  # remember to update the build number
#=================================================

build="build$V"
join="join$V.py"
base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
configdir=$builddir/configs
joinscript=$base/$join


#for mcfile in `ls $configdir | grep -E "copper|surf|bulk|reflect"`
for mcfile in `ls $configdir | grep c9- `
#for mcfile in `ls $configdir`
#for mcfile in c6-nai-surf-0.1umPb210_GRND
#for mcfile in c6-teflon-surf-1umPb210_GRND
do
    #echo -e $mcfile
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
    
    qsub -q "short" $file
    
done

#  "very_short"   #( 2 hours)
#  "short"        #(12 hours)
#  "medium"       #(48 hours)
#  "long"         #(96 hours)

