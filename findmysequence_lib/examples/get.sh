#! /bin/bash

pdbs="5me2.pdb 5irz.pdb"
maps="emd_3488.map.gz emd_8118.map.gz"

if [ -x "$(command -v wget)" ]; then
    echo "Using wget"

    for name in $maps; do
        eid=$(echo $name | awk '{print substr($1,5,4)}')
        wget ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-$eid/map/$name
    done

    for name in $pdbs; do
        wget --no-check-certificate https://files.rcsb.org/download/$name -O $name
    done

else
    echo "Using curl"

    for $name in $maps; do
        eid=$(echo $name | awk '{print substr($1,5,4)}')
        curl ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-$eid/map/$name -o $name
    done

    for $name in $pdbs; do
        curl https://files.rcsb.org/download/$name -o $name
    done
fi
gzip -d *.gz
