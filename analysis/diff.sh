#!/bin/bash
gamma=1
seedDisplacement=0
while [[ $gamma -lt 300 ]]
do
	temp=300
	while [[ $temp -lt 315 ]]
	do
		./diffusion "$temp" 0.02 1000 "$gamma" "$seedDisplacement" >> diff.dat
		seedDisplacement=$((seedDisplacement+1))
		temp=$((temp+25))
		#sleep 2
		echo "$gamma"
	done
	gamma=$((gamma+1))
done
