#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
for f in $(find . -xtype f -name "*.vcf.gz" | sort -V);do echo $f;zcat $f | grep -v "^##" | md5sum;done | paste - -
