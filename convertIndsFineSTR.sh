#!/bin/bash

# usage convertIndsFineSTR.sh firstIndIndex lastIndIndex file-prefix file-index
# converts vcf or vcf.gz to finestructure format (without header)
# to be used with vcf2fineSTR.lsf

file=$3

if [ -f $file.vcf ]; then
 for (( ind=$1+9; ind<=$2+9; ind++ ))
 do
  grep -v "^#" $file.vcf | cut -f $ind | cut -d"|" -f 1 | sed -e ':a;N;$!ba;s/\n/ /g' -e 's/ //g' >> $file.fineSTRphase.$4
  grep -v "^#" $file.vcf | cut -f $ind | cut -d"|" -f 2 |  sed -e ':a;N;$!ba;s/\n/ /g' -e 's/ //g' >> $file.fineSTRphase.$4
 done

elif [ -f $file.vcf.gz ]; then
 for (( ind=$1+9; ind<=$2+9; ind++ ))
 do
  zgrep -v "^#" $file.vcf.gz | cut -f $ind | cut -d"|" -f 1 | sed -e ':a;N;$!ba;s/\n/ /g' -e 's/ //g' >> $file.fineSTRphase.$4
  zgrep -v "^#" $file.vcf.gz | cut -f $ind | cut -d"|" -f 2 |  sed -e ':a;N;$!ba;s/\n/ /g' -e 's/ //g' >> $file.fineSTRphase.$4
 done

else
 echo "file $file.vcf[.gz] does not exist!"

fi
