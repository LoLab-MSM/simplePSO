from pysb import *


Model()
#######
V = 10.
#######
Parameter('kf1',   1./V)
Parameter('kr1',   10.)
Parameter('kcat1', 100.)
Parameter('kf2',   1./V)
Parameter('kr2',   1000.)
Parameter('kcat2', 10.)


Monomer('E', ['s'])
Monomer('S', ['e', 'type'], {'type': ['1', '2']})
Monomer('P')

# Rules
Rule('ReversibleBinding_1', E(s=None) + S(e=None, type='1') <> E(s=1) % S(e=1, type='1'), kf1, kr1)
Rule('Production_1', E(s=1) % S(e=1, type='1') >> E(s=None) + P(), kcat1)
Rule('ReversibleBinding_2', E(s=None) + S(e=None, type='2') <> E(s=1) % S(e=1, type='2'), kf2, kr2)
Rule('Production_2', E(s=1) % S(e=1, type='2') >> E(s=None) + P(), kcat2)

# Macro
# catalyze_state(E(), 's', S(), 'e', 'state', '0', '1', [kf,kr,kcat])

Observable("E_free",      E(s=None))
Observable("S1_free",     S(e=None, type='1'))
Observable("S2_free",     S(e=None, type='2'))
Observable("ES1_complex", E(s=1) % S(e=1, type='1'))
Observable("ES2_complex", E(s=1) % S(e=1, type='2'))
Observable("Product",     P())

Parameter("Etot", 1.*V)
Initial(E(s=None), Etot)

Parameter('S1_0', 10.*V)
Initial(S(e=None, type='1'), S1_0)

Parameter('S2_0', 10.*V)
Initial(S(e=None, type='2'), S2_0)
