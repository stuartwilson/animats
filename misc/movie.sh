#!/bin/bash

touch $1.mp4

ffmpeg -r 20 -i $1%05d.png -vb 5MB -vcodec mpeg4 $1.mp4 -y

open $1.mp4
