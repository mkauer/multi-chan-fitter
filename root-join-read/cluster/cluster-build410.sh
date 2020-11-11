#!/bin/bash

# 2020-08-07

#=================================================
V=410  # remember to update the build number
#=================================================

build="build$V"
join="join$V.py"
base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
configdir=$builddir/configs
joinscript=$base/$join


#for mcfile in `ls $configdir | grep -E "copper|surf|bulk|reflect"`
for mcfile in `ls $configdir | grep -E "Na24|I128"`
#for mcfile in `ls $configdir | grep c4- `
#for mcfile in `ls $configdir | grep gamma `
#for mcfile in `ls $configdir`
#for mcfile in `ls $configdir | grep data `
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
    
    #qsub -q "short" $file  ## for MC
    qsub -q "medium" $file  ## for data / and surface pb210
    #qsub -q "long" $file  ## for teflon surface expo profile
    
done

#  "very_short"   #( 2 hours)
#  "short"        #(12 hours)  good for MC
#  "medium"       #(48 hours)  good for data
#  "long"         #(96 hours)

