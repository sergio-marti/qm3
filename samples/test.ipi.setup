cat > input.xml << EOD
<simulation verbosity='high'>
  <ffsocket mode='unix' name='qm3' pbc='False'>
    <address>zundel</address>
  </ffsocket>
  <total_steps>100</total_steps>
  <output prefix='simulation'>
    <properties stride='1'> [ step, time{picosecond}, temperature{kelvin}, kinetic_cv{j/mol}, potential{j/mol} ] </properties>
    <trajectory stride='1' filename='pos'> positions{angstrom} </trajectory>
    <trajectory stride='1' filename='cen'> x_centroid{angstrom} </trajectory>
  </output>
  <system>
    <forces>
      <force forcefield='qm3'/>
    </forces>
    <initialize nbeads='4'>
      <file mode='xyz'>init.xyz</file>
      <cell mode='abc'> [1.0, 1.0, 1.0] </cell>
    </initialize>
    <ensemble>
      <temperature units='kelvin'> 300.0 </temperature>
    </ensemble>
    <motion mode='dynamics'>
      <dynamics mode='nvt'>
        <timestep units='femtosecond'> 0.25 </timestep>
        <thermostat mode='langevin'>
          <tau units='femtosecond'>100</tau>
        </thermostat>
      </dynamics>
    </motion>
  </system>
</simulation>
EOD

cat > init.xyz << EOD
7
positions{angstrom}
 O   0.54044984       -0.97485007       -0.21657970
 O   2.92484146       -0.83925223        0.15239669
 H   0.18385954       -1.25804198       -1.07875142
 H   1.74328720       -0.90927800       -0.10039502
 H   3.29706729       -1.59324218        0.64725255
 H   3.53032937       -0.60813812       -0.57607828
 H   0.04233158       -0.19485993        0.09228101
EOD
