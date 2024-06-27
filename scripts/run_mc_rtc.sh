#!/bin/bash

readonly this_dir=`cd $(dirname $0); pwd`
readonly etc_dir=${this_dir}/etc
readonly app_dir=`cd $this_dir && cd .. && cd .. && cd ..; pwd`
readonly urdf_dir="$app_dir/urdf_models"
readonly examples_dir="$app_dir/examples"


# Start mc_rtc client
MCUDPControl -s -f $etc_dir/simple.yaml &
# Get that PID
mc=$!
# Wait for the user to close AGX window
wait $agx
# Interrupt mc_rtc after AGX has been closed
kill -9 $mc
