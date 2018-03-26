from pysb import *
from pysb.macros import equilibrate

Model()

Monomer('MKK4', ['b', 'state'], {'state': ['U', 'P']})
Monomer('MKK7', ['b', 'state'], {'state': ['U', 'P']})
Monomer('JNK3', ['b', 'threo', 'tyro'], {'threo': ['U', 'P'], 'tyro': ['U', 'P']})

# All kf parameters are in units of inverse microM*s
# All kr parameters are in units of inverse s
# All kcat parameters are units of inverse s
# The forward reaction is association; the reverse is disassociation

# Because the Ordinary differential equations derived from the mass action kinetics law requires rate
# constants instead of K_D values, K_D values are going to be converted into rate parameters (k_r/k_f).
# We are going to assume that the reaction k_f is difussion limited, whereas the k_r would be allowd to vary

###### IMPORTANT INFO ABOUT JNK3


# uMKK4 with Arrestin-3, K_D = 23 microM, figure 1.A
Parameter('kf_uMKK4_Arr', 1.5e4)
Parameter('kr_uMKK4_Arr', 345000)

# pMKK4 with Arrestin-3, K_D = 347 microM, figure 1.B
Parameter('kf_pMKK4_Arr', 1.5e4)
Parameter('kr_pMKK4_Arr', 5205000)

# uMKK7 with Arrestin-3, K_D = 6.5 microM, figure 1.C
Parameter('kf_uMKK7_Arr', 1.5e4)
Parameter('kr_uMKK7_Arr', 97500) # Experimental value 97500

# pMKK7 with Arrestin-3, K_D = 13 microM, figure 1.D
Parameter('kf_pMKK7_Arr', 1.5e4)
Parameter('kr_pMKK7_Arr', 195000)

# Arrestin3-MKK4 bind to uuJNK3, K_D = 1.4 microM, figure 1.E
Parameter('kf_MKK4_ArrBinduuJNK3', 1.5e4)
Parameter('kr_MKK4_ArrBinduuJNK3', 21000)

# Arrestin3-MKK4 bind to upJNK3, K_D = 4.2 microM, figure 1.F
Parameter('kf_upJNK3Arr_MKK4', 1.5e4)
Parameter('kr_upJNK3Arr_MKK4', 63000)

# Arrestin3-MKK7 bind to puJNK3, K_D = 10.5 microM, figure 1.G
Parameter('kf_puJNK3Arr_MKK7', 1.5e4)
Parameter('kr_puJNK3Arr_MKK7', 157500)

# # upJNK3 binding arrestin, K_D = 4.2 microM, Figure 1.F
# Parameter('kf_MKK4_ArrBindupJNK3', 1.5e4)
# Parameter('kr_MKK4_ArrBindupJNK3',63000)
#
# # puJNK3 binding arrestin, K_D = 10.5 microM, Figure 1.G
# Parameter('kf_MKK4_ArrBindupJNK3', 1.5e4)
# Parameter('kr_MKK4_ArrBindupJNK3',157500)

# # ppJNK3 with Arrestin-3, K_D = 220 microM, figure 1.H
Parameter('kf_ppJNK3_Arr', 1.5e4)
Parameter('kr_ppJNK3_Arr', 3300000)

# From Pleinis et al 2017 protein expression and purification
# MKK4 activation
# MKK4: Km = 7.25 +- 1.24(microM), kcat = 8.80 +- 2.20 (inverse min) --> kcat = 0.14 (inverse s)
# MKK4: Km = (k_r + kcat) / kf, assuming difussion limited --> kr = 108749.86
Parameter('kcat_uMKK4_to_pMKK4', 0.14)

# MKK7 activation
# MKK7: Km = 3.81 +- 0.89(microM), kcat = 51.60 +- 5.80 (inverse min) --> kcat = 1.72 (inverse s)
# MKK7: Km = (k_r + kcat) / kf, assuming difussion limited --> kr = 57148.28
Parameter('kcat_uMKK7_to_pMKK7', 1.72)

# uuJNK3 binds Arrestin, K_D = 1.4 microM, figure 1.E
Parameter('kf_JNK3_Arr', 1.5e4)
Parameter('kr_JNK3_Arr', 21000)

##### FIXME: These are all guessed values

# Guessed value
# Arrestin3-MKK7 bind to uuJNK3, K_D = 1.4 microM, Figure 1.E
Parameter('kf_MKK7_ArrBinduuJNK3', 1.5e4)
Parameter('kr_MKK7_ArrBinduuJNK3', 21000)

# uJNK3 with MKK4
Parameter('kf_MKK4_uJNK3', 1.5e4)
Parameter('kr_MKK4_uJNK3', 43500)

# Guessed value
# uJNK3 with MKK7
# This is when MKK7 is bound
Parameter('kf_MKK7_uJNK3', 1.5e4)
Parameter('kr_MKK7_uJNK3', 43500)

# Equilibration between unphosphorilated MKK4 and 7, these parameters have units of inverse s
Parameter('keq_uMKK4_to_uMKK7', 10) # 7
Parameter('keq_uMKK7_to_uMKK4', 100) # 4

# # Equilibration between phosphorilated MKK4 and 7, these parameters have units of inverse s
# Parameter('kf_pMKK4_to_pMKK7', 7500)
# Parameter('kf_pMKK7_to_pMKK4', 4500)

# MKK4, Arrestin JNK3 activation
Parameter('kcat_pMKK4_ArrJNK3', 10)

# MKK7, Arrestin JNK3 activation
Parameter('kcat_pMKK7_ArrJNK3', 10)

# MKK4, JNK3 activation
Parameter('kcat_pMKK4_JNK3', 10)

# MKK7, JNK3 activation
Parameter('kcat_pMKK7_JNK3', 10)

Parameter('keq_pMKK4_to_pMKK7', 100)
Parameter('keq_pMKK7_to_pMKK4', 100)

Parameter('kf_pJNK3_MKK4complex', 1.5e4)
Parameter('kr_pJNK3_MKK4complex', 33000)

Parameter('kf_pJNK3_MKK7complex', 1.5e4)
Parameter('kr_pJNK3_MKK7complex', 33000)

# Initial conditions
Parameter('Arrestin_0', 20)
Parameter('pMKK4_0', 0.2)
Parameter('pMKK7_0', 0.2)
Parameter('uJNK3_0', 2)

Initial(MKK4(b=None, state='P'), pMKK4_0)
Initial(MKK7(b=None, state='P'), pMKK7_0)
Initial(JNK3(b=None, threo='U', tyro='U'), uJNK3_0)

# Rules

Rule('MKK4BindJNK3', MKK4(b=None, state='P') + JNK3(b=None, tyro='U') |
     MKK4(b=1, state='P') % JNK3(b=1, tyro='U'), kf_MKK4_uJNK3, kr_MKK4_uJNK3)

Rule('MKK4catJNK3', MKK4(b=1, state='P') % JNK3(b=1, tyro='U') >>
     MKK4(b=1, state='P') % JNK3(b=1, tyro='P'),kcat_pMKK4_JNK3)

# JNK3 has to get dissociated because otherwise it wouldnt be possible to have more pJNK3 than the value of MKK4 or MKK7
Rule('pJNK3_MKK4complex_diss', MKK4(b=None, state='P') + JNK3(b=None, tyro='P') |
     MKK4(b=1, state='P') % JNK3(b=1, tyro='P'), kf_pJNK3_MKK4complex, kr_pJNK3_MKK4complex)

Rule('MKK7BindJNK3', MKK7(b=None, state='P') + JNK3(b=None, threo='U') |
     MKK7(b=1, state='P') % JNK3(b=1, threo='U'), kf_MKK7_uJNK3, kr_MKK7_uJNK3)

Rule('MKK7catJNK3', MKK7(b=1, state='P') % JNK3(b=1, threo='U') >>
     MKK7(b=1, state='P') % JNK3(b=1, threo='P'), kcat_pMKK7_JNK3)

# JNK3 has to get dissociated because otherwise it wouldnt be possible to have more pJNK3 than the value of MKK4 or MKK7
Rule('pJNK3_MKK7complex_diss', MKK7(b=None, state='P') + JNK3(b=None, threo='P') |
     MKK7(b=1, state='P') % JNK3(b=1, threo='P'), kf_pJNK3_MKK7complex, kr_pJNK3_MKK7complex)


Observable('pTyr_jnk3', JNK3(b=None, tyro='P'))
Observable('pThr_jnk3', JNK3(b=None, threo='P'))
