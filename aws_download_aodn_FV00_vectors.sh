
# call as ./aws_download_aodn_FV00_vectors.sh  SYSTEM SITE YEAR STARTMONTH ENDMONTH 
# updated 20190313 to use aws instead of s3cmd to download from aodn s3 
# 20190501  modified to download entire years from AODN thredds
# 20210407 adapted FV01 yearly scipot to download FV00 monthly
# SYSTEM=  WERA or CODAR

SYSTEM=$1
SITE=$2
YEAR=$3
STARTMONTH=$4
ENDMONTH=$5

echo $SITE
echo $YEAR
echo $STARTMONTH
echo $ENDMONTH

#STARTMONTH=$3
#ENDMONTH=$4

#for month in $(seq -w $

#for i in {$STARTMONTH..$ENDMONTH}

for ((num= $STARTMONTH ; num<= $ENDMONTH; num++ ))
do
  printf -v mo "%02d" $num
echo $YEAR
echo $mo

#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 
#echo $command
echo downloading s3://imos-data/IMOS/ACORN/gridded_1h-avg-current-map_non-QC/$SITE/$YEAR/$mo
echo to
echo /home2/yasha/radar/data/$SYSTEM/RTdata_FV00/$SITE/$YEAR/$mo
mkdir -p /home2/yasha/radar/data/$SYSTEM/RTdata_FV00/$SITE/$YEAR/$mo

aws s3 sync --no-sign s3://imos-data/IMOS/ACORN/gridded_1h-avg-current-map_non-QC/$SITE/$YEAR/$mo  /home2/yasha/radar/data/$SYSTEM/RTdata_FV00/$SITE/$YEAR/$mo/ 

for day  in $(seq -w 1 31)
    do
	fprintf -v day "%02d" $dd
	echo $day
	mv /home2/yasha/radar/data/$SYSTEM/RTdata_FV00/$SITE/$YEAR/$mo/$day/*.nc /home2/yasha/radar/data/$SYSTEM/RTdata_FV00/$SITE/$YEAR/$mo/
	rmdir /home2/yasha/radar/data/$SYSTEM/RTdata_FV00/$SITE/$YEAR/$mo/$day/
done

done
#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 

#s3cmd ls s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ % this is to see the directories and what exists
# now to download it ( could do many at one time but 



#s3cmd sync s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ /home2/yasha/radar/data/WERA/radial/FRE/2018/05/