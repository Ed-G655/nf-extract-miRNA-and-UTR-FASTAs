#!/usr/bin/env bash

#Find all *out files
#find: -L option to include symlinks
find -L . \
  -type f \
  -name "*.fa.consensus" \
| sed "s#.fa.consensus#.mirmut#" \
| xargs mk
