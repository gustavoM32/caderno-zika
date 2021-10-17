#!/bin/bash

if [[ -z $1 ]]
then
	echo "Usage: ./run.sh exec"
	exit
fi

make $1

i=1

while [[ -e "in$i.txt" ]]
do
	echo "===== Test $i ====="
	echo "-- in$i.txt --"
	cat "in$i.txt"
	./$1 < "in$i.txt" > "out$i.txt"
	if [[ -e "ans$i.txt" ]]
	then
		echo "-------------"
		echo
		echo "Answer                                | Output"
		if diff -ywBW 80 --color "ans$i.txt" "out$i.txt"
		then
			echo
			echo "Passed!"
		fi
	else
		echo
		echo "-- out$i.txt --"
		cat "out$i.txt"
	fi
	echo

	i=$(($i + 1))
done