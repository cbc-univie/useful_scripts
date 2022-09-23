#!/bin/bash

# Usage: ./translate.awk yourfile.cor basefornewfile

basename=$2

# separate coordinates and header
egrep "^[ ].+[0-9].+[ ].+[0-9]" $1 > /tmp/tmpcor
egrep -v "^[ ].+[0-9].+[ ].+[0-9]" $1 | sed /^[[:space:]]*$/d > $basename.cor

# subtract minimum of each column from all values in that column
paste /tmp/tmpcor <(awk 'NR==1{a[NR]=min=$5} {a[NR]=$5; a[NR]<min?min=a[NR]:min; } END{for(i=1;i<=NR;i++){printf "%4.10f\n", a[i]-min}}' /tmp/tmpcor) <(awk 'NR==1{a[NR]=min=$6} {a[NR]=$6; a[NR]<min?min=a[NR]:min; } END{for(i=1;i<=NR;i++){printf "%4.10f\n", a[i]-min}}' /tmp/tmpcor) <(awk 'NR==1{a[NR]=min=$7} {a[NR]=$7; a[NR]<min?min=a[NR]:min; } END{for(i=1;i<=NR;i++){printf "%4.10f\n", a[i]-min}}' /tmp/tmpcor) > /tmp/newcor

# replace old coordinates by translated coordinates while preserving formatting, stolen from
# https://stackoverflow.com/questions/20835437/how-to-preserve-the-original-whitespace-between-fields-in-awk
awk 'BEGIN { FPAT = "([[:space:]]*[[:alnum:][:punct:][:digit:]]+)"; OFS = ""; } { len = length($5); $5 = sprintf("%"(len)".10f", $11); len = length($6); $6 = sprintf("%"(len)".10f", $12); $7 = sprintf("%"(len)".10f", $13); $11=$12=$13=$14=""; print $0; }' /tmp/newcor >> $basename.cor

rm /tmp/tmpcor
rm /tmp/newcor
