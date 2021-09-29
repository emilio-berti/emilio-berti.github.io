#!/bin/bash

url=https://en.wikipedia.org/wiki/$1
file=/tmp/wiki.html
usagehtml=/tmp/tmp_use.html
usagetxt=/tmp/tmp_use.txt

curl -s $url -o $file #download wikipage into temporary folder
name=$(grep -i '<title>' $file | cut -d '>' -f 2 | cut -d '<' -f 1 | cut -d '-' -f 1)
grep 'used as' $file | grep food > $usagehtml
pandoc $usagehtml -o $usagetxt #use pandoc to convert html to txt
usage=$(grep 'used as' $usagetxt | cut -d '[' -f 2 | cut -d ']' -f 1)

rm $file $usagehtml $usagetxt #remove temporary files
echo $1, $name, $usage, $url #display results