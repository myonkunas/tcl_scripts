#!/bin/bash 
catdcd -o all_.dcd control_MD29.dcd control_MD3[0-9].dcd
vmd -dispdev text -e ~/SCRIPTS/center_frames.tcl
