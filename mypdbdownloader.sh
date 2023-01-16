#!/bin/bash

# Script to download files from RCSB http file download services.
# Use the -h|--help switch to get help on usage.

# version 0.1.0     : updated from RCSB official codes, Sep 16th, 2022
# version 0.2.0     : take care of delimiter
# version 0.3.0     : add `-g' option


if ! command -v curl &> /dev/null; then
    echo "Fatal: missing command: curl"
    exit 1
fi

PROGNAME=$(basename $0)
BASE_URL='https://files.rcsb.org/download'


usage() {
    cat << EOF >&2

Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]
Usage: $PROGNAME ID,ID ID      {comma, or whitespace}

 -f <file>: the input file containing a comma or whitespace separated list of PDB ids
 -o  <dir>: the output dir, default: current dir
 -c       : download a cif[.gz] file for each PDB id
 -p       : download a pdb[.gz] file for each PDB id (default, but not available for large structures)
 -a       : download a pdb1[.gz] file (1st bioassembly) for each PDB id (not available for large structures)
 -x       : download a xml[.gz] file for each PDB id
 -s       : download a sf.cif[.gz] file for each PDB id (diffraction only)
 -m       : download a mr[.gz] file for each PDB id (NMR only)
 -r       : download a mr.str[.gz] for each PDB id (NMR only)
 -g       : whether download as the `.gz' file
EOF
    exit
}


download() {
    url="$BASE_URL/$1"
    out=$2/$1
    if [[ -f $out ]]; then
        echo "Downloaded $out, ignoring..."
    else
        echo "Downloading $url to $out"
        curl -s -f $url -o $out || echo "Failed to download $url"
    fi
}


listfile=''
outdir='.'
cif=false
pdb=false
pdb1=false
xml=false
sf=false
mr=false
mrstr=false
gz=false
ids=''
while (( $# > 0)); do
    case $1 in
        -f) listfile=$2; shift; shift; ;;
        -o) outdir=$2; shift; shift; ;;
        -c) cif=true; shift;  ;;            # care, only shift once
        -p) pdb=true; shift;  ;;
        -a) pdb1=true; shift;  ;;
        -x) xml=true; shift;  ;;
        -s) sf=true; shift;  ;;
        -m) mr=true; shift;  ;;
        -r) mrstr=true; shift;  ;;
        -h|--help) usage; exit; ;;          # exit on usage
        -g) gz=true; shift;  ;;
        *) ids="$ids $1"; shift; ;;         # shift only once!
    esac
done


# add default choice
if $cif || $pdb || $pdb1 || $xml || $sf || $mr || $mrstr; then
    :
else
    pdb=true
fi


if  [[ -f $listfile ]]; then
    contents=$(cat $listfile | tr ',' '\n' | tr ';' '\n' | tr '\n' ' ')
    for token in $contents; do
        if [[ ! $token =~ ^[0-9][a-zA-Z0-9]{3}$ ]]; then
            echo "Not a valid PDB ID: ignoring: $token"
            continue
        fi
        if $cif; then fd=${token}.cif; fi
        if $sf; then fd=${token}-sf.cif; fi
        if $pdb; then fd=${token}.pdb; fi
        if $pdb1; then fd=${token}.pdb1; fi
        if $xml; then fd=${token}.xml; fi
        if $mr; then fd=${token}.mr; fi
        if $mrstr; then fd=${token}_mr.str; fi

        if $gz; then fd=${fd}.gz; fi
        download ${fd} $outdir
    done
else
    echo "Not a file: ignoring: $listfile"
fi


ids=$(echo $ids | tr ',' '\n' | tr '\n' ' ')
for token in $ids; do
    if [[ ! $token =~ ^[0-9][a-zA-Z0-9]{3}$ ]]; then
        echo "Not a valid PDB ID: ignoring: $token"
        continue
    fi

    if $gz; then fd=${token}.pdb.gz; else fd=${token}.pdb; fi
    download ${fd} $outdir
done




