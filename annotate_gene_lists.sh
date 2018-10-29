#!/bin/bash

# -w tells grep to match whole words only (i.e. so ABC123 won't also match ABC1234).
# -F search for fixed strings (plain text) rather than regular expressions
# -f genelist.txt read search patterns from the file

GENELIST=$1
ANNOTATION=$2
OUT=$3

grep -w -F -f $GENELIST $ANNOTATION > $OUT