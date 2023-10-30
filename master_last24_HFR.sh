#!/bin/bash

# script to get recent 24 HFR nc filesplot maps, copy to webserver
#yasha hetzel 12 June 2023

unset DISPLAY
MAT="/opt/matlab_R2015a/bin/matlab -nodesktop -nosplash"
MATLAB_PATH="/home/yasha/matlab/m-files/yashatools/"

cd /home2/yasha/radar/data

./driver_get_ncrcat_last24.sh

echo "Done getting recent 24 HFR vector data...now plotting in matlab"



# now plot sites with Matlab
echo "driver_plot_recent_HFR_vectors_fn_v1_rottnest ;" | ${MAT}

# now concatinate  into pdf
echo "Done plotting 24,1,6 hr means in Matlab... Now combine 24 hr into pdf"
rm /home2/yasha/radar/data/latest24/combined_recent_HFR_24hr.pdf
convert -density 150 /home2/yasha/radar/data/latest24/0*recent*v2.png /home2/yasha/radar/data/latest24/combined_recent_HFR_24hr.pdf


#copy to website
echo "Copying recent 24 hr combined vector plot pdf to web server..."
sudo cp /home2/yasha/radar/data/latest24/combined_recent_HFR_24hr.pdf /home2/yasha/public_html/EWT/plot/yasha/radar/

echo "Done..."

################################

# Now combine and copy hourly
echo "Combine hourly into pdf"
rm /home2/yasha/radar/data/latest1hr/combined_recent_HFR_1hr.pdf
convert -density 150 /home2/yasha/radar/data/latest1hr/0*recent*v2.png /home2/yasha/radar/data/latest1hr/combined_recent_HFR_1hr.pdf


#copy to hourly to website
echo "Copying recent combined hourly vector plot pdf to web server..."
sudo cp /home2/yasha/radar/data/latest1hr/combined_recent_HFR_1hr.pdf /home2/yasha/public_html/EWT/plot/yasha/radar/

echo "Done with hourly..."

################################

# Now combine and copy 6 hour
echo "Combine 6 hour into pdf"
rm /home2/yasha/radar/data/latest1hr/combined_recent_HFR_6hr.pdf
convert -density 150 /home2/yasha/radar/data/latest6hr/0*recent*v2.png /home2/yasha/radar/data/latest1hr/combined_recent_HFR_6hr.pdf


#copy to hourly to website
echo "Copying recent combined 6 hour vector plot pdf to web server..."
sudo cp /home2/yasha/radar/data/latest1hr/combined_recent_HFR_6hr.pdf /home2/yasha/public_html/EWT/plot/yasha/radar/

echo "Done with 6 hour..."