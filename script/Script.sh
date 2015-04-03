#!/bin/bash

make
make clean
rm ./data/*

DIAM=1.0
INC=1.0
for k in {1..6}
do
	sed -i s/diameter\\t[0-9,.]*/diameter\\t$DIAM/ input_HH_CN.dat
	sed -i s/output_HH_diam_[0-9,.]*.dat/output_HH_diam_$DIAM.dat/ input_HH_CN.dat
	./main_HH_Channel_Noise_LowMemory 
	DIAM=`echo $DIAM + $INC | bc`
done

