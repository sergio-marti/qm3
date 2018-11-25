import glob
import qm3.utils.free_energy

umb = qm3.utils.free_energy.umbint( glob.glob( "pmf.???.dat" ) )
umb.setup()
umb.integrate( temperature = 300. )

wha = qm3.utils.free_energy.wham( glob.glob( "pmf.???.dat" ) )
wha.setup()
wha.integrate( temperature = 300. )

try:
    import  matplotlib.pyplot as plt
    qm3.utils.free_energy.plot_data( glob.glob( "pmf.???.dat" ) )
    plt.clf()
    plt.grid( True )
    f1, = plt.plot( umb.crd, [ i / 4.184 for i in umb.pmf ], 'g-o' )
    f1.set_label( "umbint" )
    f2, = plt.plot( wha.crd, [ i / 4.184 for i in wha.pmf ], 'b-*' )
    f2.set_label( "wham" )
    plt.legend()
    plt.tight_layout()
    plt.savefig( "pmf.pdf" )
except:
    f = open( "wham", "wt" )
    for i in range( len( wha.crd ) ):
        f.write( "%20.10lf%20.10lf\n"%( wha.crd[i], wha.pmf[i] ) )
    f.close()
    f = open( "umbint", "wt" )
    for i in range( len( umb.crd ) ):
        f.write( "%20.10lf%20.10lf\n"%( umb.crd[i], umb.pmf[i] ) )
    f.close()
