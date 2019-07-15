#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	socket


# ===================================
"""
export GauExe=`pwd`/pilgrim_client.py
export GauFchk=/usr/bin/true
./pilgrim_server.py &

pilgrim.py --gather
pilgrim.py --input

cat > pif.calcs << EOD
start_meppoint __NAME_of_TS__ gaussian
[Pilgrim_gradhess] [Pilgrim_name]
[Pilgrim_geometry]           
end_meppoint
EOD

pilgrim.py --pfn  | tee log
pilgrim.py --path | tee -a log
pilgrim.py --plot
"""
# ===================================


buf = sys.stdin.readline()
sys.stdin.readlines()
sck = socket.socket( socket.AF_UNIX, socket.SOCK_STREAM )
sck.connect( "/tmp/qm3_unix" )
sck.sendall( buf + " " * ( 1024 - len( buf ) ) )
buf = ""
while( len( buf ) < 4 ):
	buf += sck.recv( 4 )
if( buf == "done" ):
	print( "normal termination" )
sck.close()
