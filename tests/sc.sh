#!/bin/sh
TESTDIR=$(dirname $0)
julia -p 4 -L $TESTDIR/Superconvergence.jl -e '@elapsed go()'
