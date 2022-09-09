#!/bin/bash

USAGE="
Enhanced du usage, works only on current path

  -> mydu.sh                # report file and directory size (as is)
  -> mydu.sh   s|sort       # sort reports
  -> mydu.sh   d|dir        # only report directory size
  -> mydu.sh   f|file       # only report file size
"

bosort=false
bodir=false
bofile=false

while (( $# )); do
    case $1 in
        s | sort) bosort=true;;
        d | dir) bodir=true;;
        f | file) bofile=true;;
        h | -h | --help ) echo "$USAGE"; exit;;
    esac
    shift
done


## to deal with name contains whitespaces
#filetxt=()
#filesize=()
#dirtxt=()
#dirsize=()
#ls | while read i; do
#    if [[ -f "$i" ]]; then
#        filetxt+=( "$(printf '%-20s' "$i")" )
#        filesize+=( $(du -sh "$i" | awk '{print $1}') )
#        
#        for ((t=0; $t<${#filetxt[*]}; t++)); do echo "-->${filetxt[$t]}<--"; done
#        echo '========='
#        echo ">>${#filetxt[*]}<<"
#        echo '========='
#        echo ''
#        echo ''
#    fi
#
#done
## empty?? bug, cannot be dealt
#echo ${#filetxt[*]}


# to deal with name contains whitespaces
IFS=$'\n'
filetxt=()
filesize=()
dirtxt=()
dirsize=()
for i in $(ls); do
    if [[ -f "$i" ]]; then
        filetxt+=( "$(printf '%-40s' "$i")" )
        filesize+=( $(du -sh "$i" | awk '{print $1}') )
    else
        dirtxt+=( "$(printf '%-40s' "$i")" )
        dirsize+=( $(du -sh "$i" | awk '{print $1}') )
    fi
done


sortlist() {
    ndxlist=$(echo $@ | tr ' ' '\n' | sort -rh)
    reflist=''
    data=($*)
    for v in ${ndxlist[*]}; do
        for ((j=0; $j<${#data[*]}; j++)); do
            if [[ ${data[$j]} == $v ]]; then
                data[$j]="NOTHING"
                reflist="$reflist $j"
                break       # important
            fi
        done
    done
    echo "$reflist"
}


boolfile=true
booldir=true
if [[ ! $bofile && ! $bodir ]]; then
    :
elif $bofile; then
    booldir=false
elif $bodir; then
    boolfile=false
fi


if $bosort; then
    filendx=$(sortlist ${filesize[*]})
    dirndx=$(sortlist ${dirsize[*]})
    if $boolfile; then
        for i in $(echo $filendx | tr ' ' '\n'); do
            echo "file: ${filetxt[$i]}    ->  ${filesize[$i]}"
        done
    fi
    if $booldir; then
        for i in $(echo $dirndx | tr ' ' '\n'); do
            echo "dir: ${dirtxt[$i]}    ->  ${dirsize[$i]}"
        done
    fi
else
    if $boolfile; then
        for ((i=0; $i<${#filetxt[*]}; i++)); do
            echo "file: ${filetxt[$i]}    ->  ${filesize[$i]}"
        done
    fi
    if $booldir; then
        for ((i=0; $i<${#dirtxt[*]}; i++)); do
            echo "dir: ${dirtxt[$i]}    ->  ${dirsize[$i]}"
        done
    fi
fi





