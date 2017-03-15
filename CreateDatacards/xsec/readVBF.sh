#!/bin/bash


for i in $(ls | grep -e VBF)
do
  echo $i
  intlist=()

  while IFS='' read -r line || [[ -n "$line" ]]; do
    if [[ "$line" == *"accum. integral"* ]]; then
      echo "Line is "${line%%VBF*}
#      intlist+=(${line% VBF*})
    fi
  done < $i

#  echo $intlist

done
