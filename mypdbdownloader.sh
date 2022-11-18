#!/bin/bash

# Script to download files from RCSB http file download services.
# Use the -h|--help switch to get help on usage.

# version 0.1.0     : updated from RCSB official codes, Sep 16th, 2022

if ! command -v curl &> /dev/null; then
    echo "Fatal: missing command: curl"
    exit 1
fi

PROGNAME=$(basename $0)
BASE_URL='https://files.rcsb.org/download'


usage() {
    cat << EOF >&2

Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]
Usage: $PROGNAME ID,ID,ID

 -f <file>: the input file containing a comma-separated list of PDB ids
 -o  <dir>: the output dir, default: current dir
 -c       : download a cif.gz file for each PDB id
 -p       : download a pdb.gz file for each PDB id (not available for large structures)
 -a       : download a pdb1.gz file (1st bioassembly) for each PDB id (not available for large structures)
 -x       : download a xml.gz file for each PDB id
 -s       : download a sf.cif.gz file for each PDB id (diffraction only)
 -m       : download a mr.gz file for each PDB id (NMR only)
 -r       : download a mr.str.gz for each PDB id (NMR only)
EOF
    exit 1
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
        *) ids="$ids $1"; shift; ;;         # shift only once!
    esac
done

if [[ -n $listfile ]]; then
    if  [[ -f $listfile ]]; then
        contents=$(cat $listfile)
    else
        echo "Not a file: ignoring: $listfile"
    fi
fi

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"

for token in "${tokens[@]}"; do
    if [[ ! $token =~ ^[0-9][a-zA-Z]{3}$ ]]; then
        echo "Not a valid PDB ID: ignoring: $token"
        continue
    fi
    if [ "$cif" == true ]; then
        download ${token}.cif.gz $outdir
    fi
    if [ "$pdb" == true ]; then
        download ${token}.pdb.gz $outdir
    fi
    if [ "$pdb1" == true ]; then
        download ${token}.pdb1.gz $outdir
    fi
    if [ "$xml" == true ]; then
        download ${token}.xml.gz $outdir
    fi
    if [ "$sf" == true ]; then
        download ${token}-sf.cif.gz $outdir
    fi
    if [ "$mr" == true ]; then
        download ${token}.mr.gz $outdir
    fi
    if [ "$mrstr" == true ]; then
        download ${token}_mr.str.gz $outdir
    fi
done


for token in $ids; do
    if [[ ! $token =~ ^[0-9][a-zA-Z]{3}$ ]]; then
        echo "Not a valid PDB ID: ignoring: $token"
        continue
    fi
    download ${token}.pdb $outdir
done



