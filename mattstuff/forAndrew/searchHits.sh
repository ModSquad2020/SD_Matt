#!/bin/bash

###searchHits.sh
#
#
#
#Note: it will allocate 2 cores for a single job


###Usage
## sh searchHits.sh example.cm fasta.fa directoryName

# Constants
CMNAME=$(echo $1 | cut -d '.' -f1)
CMNAME=$(echo $CMNAME | cut -d "/" -f2)
CMFILE=$1
FASTAFILE=$2
DIRTOWRITE=$3

echo "Running cmsearch $CMFILE $FASTAFILE"
echo "Writing to "$DIRTOWRITE""$CMNAME".txt"
if [ ! -a $DIRTOWRITE ]
then
	mkdir $DIRTOWRITE
fi

/projects/lowelab/share/bin/x86_64/infernal1.1.2/bin/cmsearch --cpu 2 $1 $2 > $DIRTOWRITE$CMNAME.txt

