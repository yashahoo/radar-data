
# call as ./download_aodn_FV00_radials.sh  SITE YEAR STARTMONTH ENDMONTH 
SITE=$1
YEAR=$2

echo $SITE
echo $YEAR

STARTMONTH=$3
ENDMONTH=$4

#for month in $(seq -w $

#for i in {$STARTMONTH..$ENDMONTH}

for ((num= $STARTMONTH ; num<= $ENDMONTH; num++ ))
do
  printf -v mo "%02d" $num
echo $mo

#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 
#echo $command
echo downloading s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/
echo to
echo /home2/yasha/radar/data/WERA/radial/$SITE/$YEAR/$mo/
mkdir -p /home2/yasha/radar/data/WERA/radial/$SITE/$YEAR/$mo/

s3cmd sync s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/  /home2/yasha/radar/data/WERA/radial/$SITE/$YEAR/$mo/ 


done

#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 

#s3cmd ls s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ % this is to see the directories and what exists
# now to download it ( could do many at one time but 



#s3cmd sync s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ /home2/yasha/radar/data/WERA/radial/FRE/2018/05/