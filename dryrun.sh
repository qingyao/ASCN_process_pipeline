#!/usr/bin/expect

set series [lindex $argv 0]
set timeout -1

log_user 0
spawn rsync -rvnc /Users/bgprocess/aroma/hg19/processed/$series/ pgweb@130.60.23.22:/volume1/arraymapMirror/arraymap/hg19/$series/
expect "pgweb@130.60.23.22's password: "
send   "mine4genes\r"
expect "\$ "
puts $expect_out(buffer)
