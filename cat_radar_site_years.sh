# concatinate all files for given year and site
# call as ./cat_radar_site_years.sh SYSTEM SITE STARTYEAR ENDYEAR 
# this file should be run from  /Volumes/DM/data/
# yh@uwa 20231030 - updated paths for geko server 

SYSTEM=$1
SITE=$2
STARTYEAR=$3
ENDYEAR=$4

echo $SYSTEM
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

infile='IMOS*'$YEAR'*T*.nc'
outfile='IMOS_ACORN_V_'$YEAR'_'$SITE'_FV01_1-hour-avg.nc'

echo $infile
echo $outfile

#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 
#echo $command
echo downloading concatinating $SITE $YEAR 
#echo to
#echo /Volumes/DM/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR/
#mkdir -p /Volumes/DM/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR/
#aws s3 sync --no-sign s3://imos-data/IMOS/ACORN/gridded_1h-avg-current-map_QC/$SITE/$YEAR/  /Volumes/DM/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR/ 

#find /Volumes/DM/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR -name 'IMOS*2007*T*.nc' -type f -print | sort | ncrcat -h IMOS_ACORN_V_$YEAR_$SITE_FV01_1-hour-avg.nc
find /Volumes/DM/data/$SYSTEM/RTdata_FV01/$SITE/$YEAR -name $infile -type f -print | sort | ncrcat -h $outfile

done

#command=s3://imos-data/IMOS/ACORN/radial/$SITE/$YEAR/$mo/ 

#s3cmd ls s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ % this is to see the directories and what exists
# now to download it ( could do many at one time but 



#s3cmd sync s3://imos-data/IMOS/ACORN/radial/FRE/2018/05/ /Volumes/DM/data/WERA/radial/FRE/2018/05/