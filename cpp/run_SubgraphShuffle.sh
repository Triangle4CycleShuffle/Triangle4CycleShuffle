#!/bin/bash -x

if [ $# -ne 1 ]; then
    echo "USAGE: run_SubgraphShuffle.sh [Dataset]"
    exit 1
fi

./SubgraphShuffle ../data/${1}/edges.csv 20000 1-8 -1 20-1 2n
./SubgraphShuffle ../data/${1}/edges.csv 20000 1-8 -1 20-1 3n
./SubgraphShuffle ../data/${1}/edges.csv 20000 1-8 -1 20-1 4n
./SubgraphShuffle ../data/${1}/edges.csv 20000 1-8 -1 20-1 5
./SubgraphShuffle ../data/${1}/edges.csv 20000 1-8 -1 20-1 6-0.1
./SubgraphShuffle ../data/${1}/edges.csv 20000 1-8 -1 20-1 8n
./SubgraphShuffle ../data/${1}/edges.csv 20000 1-8 -1 20-1 9
