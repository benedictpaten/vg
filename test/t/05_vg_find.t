#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for vg

plan tests 5

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
is $? 0 "construction"

vg index -s x.vg
is $? 0 "indexing nodes and edges of graph"

vg index -k 11 x.vg
is $? 0 "indexing 11mers"

node_matches=$(vg find -k TAAGGTTTGAA -c 0 x.vg | vg view -g - | grep "^S" | cut -f 2 | grep '1$\|2$\|9$\|5$\|6$\|8$' | wc -l)
is $node_matches 6 "all expected nodes found via kmer find"

edge_matches=$(vg find -k TAAGGTTTGAA -c 0 x.vg | vg view -g - | grep "^L" | cut -f 2 | grep '1$\|2$\|8$\|5$\|6$' | wc -l)
is $edge_matches 5 "all expected edges found via kmer find"

rm -rf x.vg.index
rm -f x.vg
