#!/bin/bash

mkdir salmon_sf
for i in *_quant;
do
  	p=$(echo $i| cut -d'_' -f 1,2);
  	cd ${p}_quant;
  	cp quant.sf ${p}_quant.sf;
  	mv ${p}_quant.sf ../salmon_sf;
  	cd ../
done
