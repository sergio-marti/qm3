import qm3.utils.free_energy

dene = [ 6.372, 5.995, 5.893, 4.233, 5.740, 5.230, 5.184, 4.128, 5.336, 5.004, 3.907, 5.164,
    6.588, 6.661, 5.443, 6.704, 5.136, 2.071, 3.673, 4.948, 6.117, 5.776, 5.247, 7.117, 4.269,
    6.832, 6.346, 6.158, 3.961, 3.197, 4.638, 5.314, 5.121, 6.549, 6.110, 5.151, 4.629, 4.772,
    5.990, 5.979, 6.092, 5.058, 5.820, 5.658, 5.913, 6.007, 5.404, 4.699, 5.271, 6.511, 7.027 ]

for dat in qm3.utils.free_energy.fep_integrate( dene ):
    print( "Samples: ", dat["Samples"] )
    print( "dF:      ", dat["dF"]/4.184, "+/-", dat["Error"]/4.184, "_kcal/mol" )
    print( "Samp.rat:", dat["Sampling Ratio"] )
    print( "Autocorr:", dat["Autocorrelation"] )

for dat in qm3.utils.free_energy.fep_integrate( dene, clusters = 2, tries = 10 ):
    print( 80 * "-" )
    print( "Samples: ", dat["Samples"] )
    print( "dF:      ", dat["dF"]/4.184, "+/-", dat["Error"]/4.184, "_kcal/mol" )
    print( "Samp.rat:", dat["Sampling Ratio"] )
    print( "Autocorr:", dat["Autocorrelation"] )
