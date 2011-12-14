#!/bin/sh

###########################################
#
# plot_gauss.sh
#
# Based on SDDSTOOLS/plot_gauss.sh
#
# makes gauss fits

# -----------------------------------------
#
# Usage
# =====
#
# plot_gauss.sh data.sdds
#
# -----------------------------------------
#
# Arguments
# =========
#
# data.sdds : SDDS file from which the data
#             is to be extracted. Expected
#             to be output from Elegant.
#             If the file's extension is
#             either '.sdds', or '.out',
#             it may be omitted.
###########################################

# usage
if ! ( [ $1 ] && [ $2 ] && [ $3 ] ); then
    echo "usage:"
    echo "plot_gauss.sh data.sdds"
    exit 1
fi

# define variables from arguments
data_base=${1%.*}

# get full data file name
if [ "$1" == "$data_base" ];then
    if [ -e $data_base.out ]; then
	data=$data_base.out
    fi
    if [ -e $data_base.sdds ]; then
	data=$data_base.sdds
    fi
else
    data=$1
fi

# unit conversion and definition of columns
sddsprocess $data $data_base.gauss -noWarnings \
 "-convertUnits=column,x,mm,m,1.e3" \
 "-convertUnits=column,xp,mrad,,1.e3" \
 "-convertUnits=column,y,mm,m,1.e3" \
 "-convertUnits=column,yp,mrad,,1.e3" \
 "-process=t,average,tbar" \
 "-define=column,z,t tbar - 2.99792458E11 *,units=mm" \
 "-define=column,d,p pCentral - pCentral / 100 *,units=%"
# calculate rms values

sddsprocess $data_base.gauss -noWarnings \
 "-process=x,standardDeviation,xrms" \
 "-process=xp,standardDeviation,xprms" \
 "-process=y,standardDeviation,yrms" \
 "-process=yp,standardDeviation,yprms" \
 "-process=z,standardDeviation,zrms" \
 "-process=d,standardDeviation,drms" \
 "-process=t,standardDeviation,trms" \
 "-process=p,standardDeviation,prms"


# create histograms
sddshist $data_base.gauss $data_base.xhis  -data=x  -bin=100
sddshist $data_base.gauss $data_base.xphis -data=xp -bin=100
sddshist $data_base.gauss $data_base.yhis  -data=y  -bin=100
sddshist $data_base.gauss $data_base.yphis -data=yp -bin=100
sddshist $data_base.gauss $data_base.zhis  -data=z  -bin=100
sddshist $data_base.gauss $data_base.dhis  -data=d  -bin=100
sddshist $data_base.gauss $data_base.this  -data=t  -bin=100
sddshist $data_base.gauss $data_base.phis  -data=p  -bin=100

# perform gaussian fits
sddsgfit $data_base.xhis  xhis.gfit  -column=x,frequency
sddsgfit $data_base.xphis xphis.gfit -column=xp,frequency
sddsgfit $data_base.yhis  yhis.gfit  -column=y,frequency
sddsgfit $data_base.yphis yphis.gfit -column=yp,frequency
sddsgfit $data_base.zhis  zhis.gfit  -column=z,frequency
sddsgfit $data_base.dhis  dhis.gfit  -column=d,frequency
sddsgfit $data_base.this  this.gfit  -column=t,frequency
sddsgfit $data_base.phis  phis.gfit  -column=p,frequency

# clean up
rm $data_base.gauss
rm $data_base.gauss~
rm $data_base.*his
rm $data_base.*his~
rm *his.gfit~
#rm *his.gfit
