#!/bin/bash

###doSearch.sh
#Purpose is to CMSEARCH through a fasta file
#against a directory of CMs.

#Give this a directory of CMs
#The fasta file of sequences to search, and
#The new directory to write to.

###Usage
# sh doSearch.sh dirWithCMs fastaFile.fa dirToWrite

###Constants
path=$1
FASTAFILE=$2
DIRTOWRITE=$3

MAXJOBS=17

for f in ${path}*.cm
do
	echo $f
	sh searchHits.sh $f $FASTAFILE $DIRTOWRITE &

	# https://stackoverflow.com/questions/8734998/how-can-i-throttle-the-number-of-commands-run-inside-a-bash-script?noredirect=1&lq=1
	while true
	do
		NUMJOBS=`ps --ppid $$ -o pid= | wc -l`
		test $NUMJOBS -lt $MAXJOBS && break
	done
done
echo "Done running doSearch.sh"


### Cleaning the directory of weird .xxx artifacts
rm ${DIRTOWRITE}*.xxx
echo "Done cleaning .xxx artifacts"
