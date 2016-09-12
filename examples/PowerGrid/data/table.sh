#listing the summary of the power grid data
echo "instance name \\ description \\ buses \\ branches \\ generators"
for f in `ls case*.m`
do
    casename=`basename $f .m`
    CASENAME=`echo $casename | awk '{print toupper($0)}'`
    desc=`grep $CASENAME $f | cut -d' ' -f2-`
    nbus=`cat ${casename}.bus | wc -l`
    nbranch=`cat ${casename}.branch | wc -l`
    ngen=`cat ${casename}.gen | wc -l`
    echo "$casename\\$desc\\$nbus\\$nbranch\\$ngen"
done
