#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import  sys
import  os
import  time
import  glob
import  random
import  socket
import  struct
import  time

xbx = [ 0, 29 ]
ybx = [ 0, 43 ]
chg = [ (-1,0), (+1,0), (0,-1), (0,+1) ]

random.seed()
cur = []
for ff in glob.glob( "pes.??.??" ):
    tmp = ff.split( "." )
    cur.append( ( int( tmp[1] ), int( tmp[2] ) ) )
end = struct.pack( "i", -1 ) + struct.pack( "i", -1 ) + struct.pack( "i", -1 ) + struct.pack( "i", -1 )
run = []
ssk = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
flg = True
while( flg ):
    try:
        ssk.bind( ( "127.0.0.1", 6969 ) )
        flg = False
    except:
        time.sleep( 1 )
print( "+ ready!" )
while( True ):
    ssk.listen( 100 )
    chd, adr = ssk.accept()
    # -----------------------------------------------------------
    for i in range( len( run ) -1, -1, -1 ):
        tx, ty = run[i]
        if( os.path.isfile( "pes.%02d.%02d"%( tx, ty ) ) ):
            cur.append( ( tx, ty ) )
            del run[i]
    print( "-- updated:", len( cur ), len( run ) )
    # -----------------------------------------------------------
    lst = []
    for ox, oy in cur:
        random.shuffle( chg )
        for dx, dy in chg:
            if( ox+dx >= xbx[0] and ox+dx <= xbx[1] and
                oy+dy >= ybx[0] and oy+dy <= ybx[1] and
                not ( ox+dx, oy+dy ) in cur and not ( ox+dx, oy+dy ) in run ):
                lst.append( ( ox, oy, ox+dx, oy+dy ) )
    if( len( lst ) > 0 ):
        ox,oy,cx,cy = random.choice( lst )
        run.append( ( cx, cy ) )
        chd.send( struct.pack( "i", ox ) + struct.pack( "i", oy ) + struct.pack( "i", cx ) + struct.pack( "i", cy ) )
        print( ox, oy, ">>", cx, cy )
    else:
        chd.send( end )
    chd.close()
    # -----------------------------------------------------------
