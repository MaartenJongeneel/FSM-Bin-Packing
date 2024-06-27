#!/bin/bash

readonly this_dir=`cd $(dirname $0); pwd`
rosrun rviz rviz -d $this_dir/etc/display.rviz
