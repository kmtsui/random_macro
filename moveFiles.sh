#!/bin/sh

for (( i=0; i < 700; ++i ))
do
    newIdx=$((i+300))
    cp BnL/unbinned/diffuser4_400nm_toy${i}.root BnL/test/diffuser4_400nm_toy${newIdx}.root
done