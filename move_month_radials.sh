# call as ./move_month_radials.sh  STARTD ENDD 
# this is redundant with simone m-file but more effiecient if edit a bit more
# this is copied to processing directory for syncing with github but will not work if not in data directory
# in fact i think this file is not in use and is broken so don't think it is relevent in present workflow
# yasha 2020-02-25

STARTD=$1
ENDD=$2
MON=05d
mkdir -p $MON

#for((day=$STARTD  ; day<$ENDD+1 ; day++))
#for((dd=$STARTD  ; dd<$ENDD+1 ; dd++))
for day  in $(seq -w $STARTD $ENDD)
#fprintf -v day "%02d" $dd

#echo $day

	do


#mv $day/*.nc ./$MON/
rmdir $day

echo done copying $day

done