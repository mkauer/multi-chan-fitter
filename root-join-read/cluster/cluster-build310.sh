#!/bin/bash

# 2019-04-25

#=================================================
V=310  # remember to update the build number
#=================================================

build="build$V"
join="join$V.py"
base=/home/mkauer/mc-fitting/root-join-read
clustdir=$base/cluster
builddir=$base/$build
configdir=$builddir/configs
joinscript=$base/$join


#for mcfile in `ls $configdir | grep -E "SET1-V00-04-12|SET2-V00-04-12"`
#for mcfile in `ls $configdir | grep data`
#for mcfile in `ls $configdir`
for mcfile in `ls $configdir | grep c7-data-1616-V00-04-12`	      
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

