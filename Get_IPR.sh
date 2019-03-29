#!/bin/bash
#1.1 0.75 0.85 1.3 0.05 0.025 0.075
for dis in 1.0 #0.0001 0.0005 0.001 0.002 0.005 0.01 0.025 0.05 0.075 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.85 0.9 1.0 1.1 1.2 1.3 1.4 2.0 3.0 4.0 #5.0 #0.4 0.5 0.6 0.65 0.7 0.8 0.9 1.0 1.2 1.4 2.0 3.0 4.0 5.0
do printf "${dis} "
std_dev=0.0
mean=0.0

for conf in {1..20}
do val=$(awk 'NR == 30 {print $7}' Runs/J_1.15/Dis_${dis}/unbiased_OP_seed_101/Disorder_seed_${conf}/out_run.txt)
printf "${val} "
mean=$(echo "$mean + $val" | bc)
done
mean=$(echo "${mean}*(0.05)" | bc)
printf "${mean} "


for conf in {1..20}
do val=$(awk 'NR == 30 {print $7}' Runs/J_1.15/Dis_${dis}/unbiased_OP_seed_101/Disorder_seed_${conf}/out_run.txt)
std_dev=$(echo "${std_dev} + ((${val} - ${mean})*(${val} - ${mean}))" | bc)
done
std_dev=$(echo "${std_dev}*(0.05)" | bc)
std_dev=$(echo "sqrt(${std_dev})" | bc)


printf "${std_dev} "
echo ""
done
