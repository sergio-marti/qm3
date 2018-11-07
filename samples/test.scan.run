import os
import time
import threading


LK = threading.Lock()
Ni = 20 + 1
Nj = 20 + 1
Np = 4
FF = []
CC = []
for i in range( Ni ):
    for j in range( Nj ):
        CC.append( ( i, j ) )
        FF.append( not os.path.isfile( "pes.%d.%d"%( i, j ) ) )
NN = len( FF )


def worker( num ):
    global    LK, CC, FF, NN
    fd = open( "run.%d.log"%( num ), "wt" )
    while( sum( FF ) > 0 ):
        LK.acquire()
        w  = 0
        Wi = None
        Wj = None
        while( w < NN and ( Wi == None or Wj == None ) ):
            if( FF[w] and os.path.isfile( "pes.%d.%d"%( CC[w][0], CC[w][1] - 1 ) ) ):
                Wi = CC[w][0]; Oi = CC[w][0]
                Wj = CC[w][1]; Oj = CC[w][1] - 1
            elif( FF[w] and os.path.isfile( "pes.%d.%d"%( CC[w][0] - 1, CC[w][1] ) ) ):
                Wi = CC[w][0]; Oi = CC[w][0] - 1
                Wj = CC[w][1]; Oj = CC[w][1]
            elif( FF[w] and os.path.isfile( "pes.%d.%d"%( CC[w][0] - 1, CC[w][1] - 1 ) ) ):
                Wi = CC[w][0]; Oi = CC[w][0] - 1
                Wj = CC[w][1]; Oj = CC[w][1] - 1
            else:
                w += 1
        if( Wi != None and Wj != None ):
            FF[w] = False
        LK.release()
        if( Wi != None and Wj != None ):
            fd.write( "%d.%d >> %d.%d"%( Oi, Oj, Wi, Wj ) ); fd.flush()
            os.system( "python3 test.scan %d %d %d %d"%( Oi, Oj, Wi, Wj ) )
            fd.write( " >> done!\n" ); fd.flush()
        time.sleep( 1 )


pid = []
for i in range( Np ):
    pid.append( threading.Thread( target = worker, args = ( i, ) ) )
    pid[-1].start()
    time.sleep( 1 )
for i in range( Np ):
    pid[i].join()
