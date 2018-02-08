#!/bin/bash

# usage convertIndsFineSTR.sh firstIndIndex lastIndIndex file-prefix file-index

file=$3


for (( ind=$1+9; ind<=$2+9; ind++ ))
do
 grep -v "^#" $file.vcf | cut -f $ind | cut -d"|" -f 1 | sed -e ':a;N;$!ba;s/\n/ /g' -e 's/ //g' >> $file.fineSTRphase.$4
 grep -v "^#" $file.vcf | cut -f $ind | cut -d"|" -f 2 |  sed -e ':a;N;$!ba;s/\n/ /g' -e 's/ //g' >> $file.fineSTRphase.$4
done

