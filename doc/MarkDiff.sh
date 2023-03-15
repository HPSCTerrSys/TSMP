#!/bin/bash

origDir=$1
newDir=$2

cd ${newDir}
for file in $(find ./ -type f -name "*.html")
do
  echo "Found file: ${file}"
  if [[ -f "${origDir}/${file}" ]]; then
    diff --new-line-format='<b style="color:#b96922;">%L</b>' \
      --old-line-format='' --unchanged-line-format='%L' \
      ${origDir}/${file} ${file} > ${file}_tmp
    mv "${file}_tmp" ${file}
  else
    echo "DEBUG: Skip marking diff of ${file} as not present in ${origDir}"
  fi
done

