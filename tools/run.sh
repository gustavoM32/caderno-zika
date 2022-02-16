#!/bin/bash

if [[ -z $1 ]]
then
	echo "Usage: ./run.sh exec"
	exit
fi

if ! make $1
then
	exit
fi

for i in {1..9}
do
	if [[ -e "in$i.txt" ]]
	then
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
				echo -e "\e[1m\e[32mPassed!\e[39m\e[0m"
			else
				echo
				echo -e "\e[1m\e[31mFailed\e[39m\e[0m"
			fi
		else
			echo
			echo "-- out$i.txt --"
			cat "out$i.txt"
			echo
			echo -e "\e[1m\e[33mManual\e[39m\e[0m"
		fi
		echo
	fi
done