
# call as ./aws_download_aodn_FV01_vectors.sh  SYSTEM SITE STARTYEAR ENDYEAR 
# updated 20190313 to use aws instead of s3cmd to download from aodn s3 
# 20190501  modified to download entire years from AODN thredds
# SYSTEM=  WERA or CODAR

SYSTEM=$1
SITE=$2
STARTYEAR=$3
ENDYEAR=$4

echo $SITE
echo $YEAR

#STARTMONTH=$3
#ENDMONTH=$4

#for month in $(seq -w $

#for i in {$STARTMONTH..$ENDMONTH}

for ((num= $STARTYEAR ; num<= $ENDYEAR; num++ ))
do
  printf -v YEAR "%04d" $num
echo $YEAR

#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 
#echo $command
echo downloading s3://imos-data/IMOS/ACORN/gridded_1h-avg-current-map_QC/$SITE/$YEAR/
echo to
echo /home2/yasha/radar/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR/
mkdir -p /home2/yasha/radar/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR/

aws s3 sync --no-sign s3://imos-data/IMOS/ACORN/gridded_1h-avg-current-map_QC/$SITE/$YEAR/  /home2/yasha/radar/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR/ 


done

#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 

#s3cmd ls s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ % this is to see the directories and what exists
# now to download it ( could do many at one time but 



#s3cmd sync s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ /home2/yasha/radar/data/WERA/radial/FRE/2018/05/