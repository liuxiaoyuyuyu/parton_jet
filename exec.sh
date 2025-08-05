#!/bin/bash

#root "evolution_jet_bins.cc($1)"
root -l -b -q "parton_v2.C(\"$1\")" 