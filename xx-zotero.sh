#!/bin/bash

# version 0.1.0 : split/combine big file


usage="Usage:
    >>> xx-zotero split FILE     (results will be FILE-n, where n is integers)
    >>> xx-zotero combine FILEs  (every FILE needs end with increasing number)
"

if (( $# < 2 )); then
    echo "Fatal: wrong input"
    echo "$usage"
    exit
elif [[ $1 == split ]]; then
    if (( $# != 2 )); then echo "Fatal: split: wrong input"; echo "$usage"; exit; fi
    if [[ ! -f $2 ]]; then echo "Fatal: split: not a file: $2"; exit; fi
elif [[ $1 == combine ]]; then
    shift
    base=${1%--*}
    for f in $*; do
        if [[ ! -f $f ]]; then echo "Fatal: combine: not a file: $2"; exit; fi
        name=${f%--*}
        if [[ $name != $base ]]; then echo "Fatal: combine: file not consistent: $f"; exit; fi
        num=${f##*--}
        if [[ ! $num =~ ^[0123456789]+$ ]]; then
            echo "Fatal: combine: wrong file format: $f"
            exit
        fi
    done
else
    echo "Fatal: not supported"
    echo "$usage"
    exit
fi

folder=x-zotero-split
mkdir -p ~/$folder

if [[ $1 == split ]]; then
    echo "Note: results will be saved in folder: ~/$folder"
    file=$(realpath $2)
    base=$(basename $2)
    cd ~/$folder/
    split -d -b 2G $file $base--
    exit
fi

cd ~/$folder
# aleady `shift`ed
if [[ -f $base ]]; then
    echo "Warning: file already exist: $base"
    exit
fi
echo "Note: result: ~/$folder/$base"
cat $* > $base

