from pysb import *
from pysb.macros import equilibrate

Model()

Monomer('Arrestin', ['b1', 'b2', 'b3'])
Monomer('ASK1', ['b', 'state'], {'state': ['U', 'P']})
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


# uMKK4 with Arrestin-3, K_D = 23 microM, figure 1.A
Parameter('kf_uMKK4_Arr', 1.5e4)
Parameter('kr_uMKK4_Arr', 345000)

# # pMKK4 with Arrestin-3, K_D = 347 microM, figure 1.B
# Parameter('kf_pMKK4_Arr', 1.5e4)
# Parameter('kr_pMKK4_Arr', 5205000)

# uMKK7 with Arrestin-3, K_D = 6.5 microM, figure 1.C
Parameter('kf_uMKK7_Arr', 1.5e4)
Parameter('kr_uMKK7_Arr', 97500)

# # pMKK7 with Arrestin-3, K_D = 13 microM, figure 1.D
# Parameter('kf_pMKK7_Arr', 1.5e4)
# Parameter('kr_pMKK7_Arr', 195000)

# uJNK3 with Arrestin-3, K_D = 1.4 microM, figure 1.E
Parameter('kf_uJNK3_Arr', 1.5e4)
Parameter('kr_uJNK3_Arr', 21000)

# # pJNK3 with Arrestin-3, K_D = 220 microM, figure 1.F
# Parameter('kf_pJNK3_Arr', 1.5e4)
# Parameter('kr_pJNK3_Arr', 3300000)

# From Pleinis et al 2017 protein expression and purification
# MKK4 activation
# MKK4: Km = 7.25 +- 1.24(microM), kcat = 8.80 +- 2.20 (inverse min) --> kcat = 0.14 (inverse s)
# MKK4: Km = (k_r + kcat) / kf, assuming difussion limited --> kr = 108749.86
Parameter('kcat_uMKK4_to_pMKK4', 0.14)

# MKK7 activation
# MKK7: Km = 3.81 +- 0.89(microM), kcat = 51.60 +- 5.80 (inverse min) --> kcat = 1.72 (inverse s)
# MKK7: Km = (k_r + kcat) / kf, assuming difussion limited --> kr = 57148.28
Parameter('kcat_uMKK7_to_pMKK7', 1.72)


# FIXME: These are all guessed values

# Ask1 bind to Arrestin, no KD data
Parameter('kf_Ask1_Arr', 1.5e4)
Parameter('kr_Ask1_Arr', 300000)
Parameter('kcat_Ask1_Activation', 4)

# Equilibration between unphosphorilated MKK4 and 7, these parameters have units of inverse s
Parameter('kf_uMKK4_to_uMKK7', 7)
Parameter('kf_uMKK7_to_uMKK4', 4)

# # Equilibration between phosphorilated MKK4 and 7, these parameters have units of inverse s
# Parameter('kf_pMKK4_to_pMKK7', 7500)
# Parameter('kf_pMKK7_to_pMKK4', 4500)

# MKK4, JNK3 activation
Parameter('kcat_pMKK4_JNK3', 5)

# MKK7, JNK3 activation
Parameter('kcat_pMKK7_JNK3', 7)

# # MKK4 ASK1 release
# Parameter('kf_MKK4_Ask1_release', 5000)
#
# # MKK7 ASK1 release
# Parameter('kf_MKK7_Ask1_release', 7000)

# Initial conditions
Parameter('Ask1_0', 0.05)
Parameter('Arrestin_0', 0.05)
Parameter('uMKK4_0', 0.05)
Parameter('uMKK7_0', 0.05)
Parameter('uJNK3_0', 0.05)

Initial(ASK1(b=None, state='U'), Ask1_0)
Initial(Arrestin(b1=None, b2=None, b3=None), Arrestin_0)
Initial(MKK4(b=None, state='U'), uMKK4_0)
Initial(MKK7(b=None, state='U'), uMKK7_0)
Initial(JNK3(b=None, threo='U', tyro='U'), uJNK3_0)

# Rules

Rule('Ask1bindArrestin', ASK1(b=None, state='U') + Arrestin(b1=None, b2=None, b3=None) |
     ASK1(b=1, state='U') % Arrestin(b1=1, b2=None, b3=None), kf_Ask1_Arr, kr_Ask1_Arr)

Rule('Ask1Activation', ASK1(b=1, state='U') % Arrestin(b1=1, b2=None, b3=None) >>
     ASK1(b=1, state='P') % Arrestin(b1=1, b2=None, b3=None), kcat_Ask1_Activation)

Rule('Ask1ArrBindMkk4', ASK1(b=1, state='P') % Arrestin(b1=1, b2=None, b3=None) + MKK4(b=None, state='U') |
     ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='U'), kf_uMKK4_Arr, kr_uMKK4_Arr)

Rule('Ask1ArrMkk4Activation', ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='U') >>
            ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='P'), kcat_uMKK4_to_pMKK4)

Rule('Ask1ArrBindMkk7', ASK1(b=1, state='P') % Arrestin(b1=1, b2=None, b3=None) + MKK7(b=None, state='U') |
     ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='U'), kf_uMKK7_Arr, kr_uMKK7_Arr)

Rule('Ask1ArrMkk7Activation', ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='U') >>
            ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='P'), kcat_uMKK7_to_pMKK7)
# EquilibrateUMKK4And7
equilibrate(ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='U'),
            ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='U'), [kf_uMKK4_to_uMKK7, kf_uMKK7_to_uMKK4])

# EquilibratePMKK4and7
# equilibrate(ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='P'),
#             ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='P'), [kf_pMKK4_to_pMKK7, kf_pMKK7_to_pMKK4])

# The other hipothesis is that Ask1 gets released
# Rule('MKK4Ask1Release', ASK1(b=1) % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='P') >>
#       Arrestin(b1=None, b2=2, b3=None) % MKK4(b=2, state='P') + ASK1(b=None), kf_MKK4_Ask1_release)
# Rule('MKK7Ask1Release', ASK1(b=1) % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='P') >>
#       Arrestin(b1=None, b2=2, b3=None) % MKK7(b=2, state='P') + ASK1(b=None), kf_MKK7_Ask1_release)

Rule('MKK4BindJNK3', ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='P') + JNK3(b=None, tyro='U') |
     ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, tyro='U'), kf_uJNK3_Arr, kr_uJNK3_Arr)

Rule('MKK4catJNK3', ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, tyro='U') >>
     ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK4(b=2, state='P') + JNK3(b=None, tyro='P'), kcat_pMKK4_JNK3)

Rule('MKK7BindJNK3', ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='P') + JNK3(b=None, threo='U') |
     ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='U'), kf_uJNK3_Arr, kr_uJNK3_Arr)

Rule('MKK7catJNK3', ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='U') >>
     ASK1(b=1, state='P') % Arrestin(b1=1, b2=2, b3=None) % MKK7(b=2, state='P') + JNK3(b=None, threo='P'), kcat_pMKK7_JNK3)

Observable('mkk4_pjnk3', JNK3(b=None, threo='U', tyro='P'))
Observable('mkk7_pjnk3', JNK3(b=None, threo='P', tyro='U'))



