#!/bin/sh
#
# This is a bash file to perform QC on SIGMA-AWS data
# Motoshi Nishimura @ National Institute of Polar Research
#
###########################################################################

cyear=12-20

site_flag=SIGMA-A
#site_flag=SIGMA-B

#####
echo '===================================================='
start_time=$(date +'%s.%3N')
echo 'Program start time:'
echo "start time: $(date -d "@${start_time}" +'%Y-%m-%d %H:%M:%S.%3N (%:z)')"
echo ''


if [ ${site_flag} = SIGMA-A ] ; then
  cp ../input/SIGMA-A_${cyear}_met_Lv1.2.csv ./${cyear}.csv
fi

if [ ${site_flag} = SIGMA-B ] ; then  
  cp ../input/SIGMA-B_${cyear}_met_Lv1.2.csv ./${cyear}.csv
fi

input="${cyear}-${site_flag}"

sed -i -e "s/\//,/g" ${cyear}.csv
sed -i -e "s/ /,/g" ${cyear}.csv
sed -i -e "s/,,/,/g" ${cyear}.csv
sed -i -e "s/,,,/,/g" ${cyear}.csv
sed -i -e "s/,,/,/g" ${cyear}.csv

echo "${input}" > datalist.txt
./make_Level1.3data < datalist.txt

rm datalist.txt
rm ./${cyear}.csv

mv ./*.csv ../output

end_time=$(date +'%s.%3N')
echo 'Program end time:'
echo "end time: $(date -d "@${end_time}" +'%Y-%m-%d %H:%M:%S.%3N (%:z)')"
