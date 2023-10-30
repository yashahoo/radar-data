# call as ./copy2queue.sh  SITE YEAR STARTMONTH ENDMONTH 
# updated 20190313 to use aws instead of s3cmd to download from aodn s3 

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
#echo $mo


echo /data1/Yasha/AODN_uploads/RTdata_FV01/radial/$SITE/$YEAR$mo/*.nc /data1/Yasha/queued/


echo done copying $mo

done