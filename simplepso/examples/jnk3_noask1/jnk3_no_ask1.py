from pysb import *
from pysb.macros import equilibrate

Model()

# Declaring the monomers of the model
Monomer('Arrestin', ['b1', 'b2', 'b3'])
Monomer('MKK4', ['b', 'state'], {'state': ['U', 'P']})
Monomer('MKK7', ['b', 'state'], {'state': ['U', 'P']})
Monomer('JNK3', ['b', 'threo', 'tyro'], {'threo': ['U', 'P'], 'tyro': ['U', 'P']})

# All kf parameters are in units of inverse microM*s
# All kr parameters are in units of inverse s
# All kcat parameters are units of inverse s
# The forward reaction is association; the reverse is disassociation

# Because the Ordinary differential equations derived from the mass action kinetics law requires rate
# constants instead of K_D values, K_D values are going to be converted into rate parameters (k_r/k_f).
# We are going to assume that the reaction k_f is difussion limited, whereas the k_r would be allowed to vary

###### IMPORTANT INFO ABOUT JNK3

# uMKK4 with Arrestin-3, K_D = 23 microM, figure 1.A
Parameter('kf_uMKK4_Arr', 2) # 1.5e4
Parameter('kr_uMKK4_Arr', 46)

# pMKK4 with Arrestin-3, K_D = 347 microM, figure 1.B
Parameter('kf_pMKK4_Arr', 2)
Parameter('kr_pMKK4_Arr', 80)

# uMKK7 with Arrestin-3, K_D = 6.5 microM, figure 1.C
Parameter('kf_uMKK7_Arr', 2)
Parameter('kr_uMKK7_Arr', 13) # Experimental value 97500

# pMKK7 with Arrestin-3, K_D = 13 microM, figure 1.D
Parameter('kf_pMKK7_Arr', 2)
Parameter('kr_pMKK7_Arr', 26)

# Arrestin3-MKK4 bind to uuJNK3, K_D = 1.4 microM, figure 1.E
Parameter('kf_MKK4_ArrBinduuJNK3', 2)
Parameter('kr_MKK4_ArrBinduuJNK3', 2.8)

# Arrestin3-MKK4 bind to upJNK3, K_D = 4.2 microM, figure 1.F
Parameter('kf_upJNK3Arr_MKK4', 2)
Parameter('kr_upJNK3Arr_MKK4', 8.4)

# Arrestin3-MKK7 bind to puJNK3, K_D = 10.5 microM, figure 1.G
Parameter('kf_puJNK3Arr_MKK7', 2)
Parameter('kr_puJNK3Arr_MKK7', 21)

# # ppJNK3 with Arrestin-3, K_D = 220 microM, figure 1.H
Parameter('kf_ppJNK3_Arr', 2)
Parameter('kr_ppJNK3_Arr', 440)

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
Parameter('kf_JNK3_Arr', 2)
Parameter('kr_JNK3_Arr', 2.8)

##### These are the parameters that are going to be calibrated

Parameter('kf_MKK4BindArr_JNK3', 2)
Parameter('kr_MKK4BindArr_JNK3', 2)

Parameter('kf_MKK7BindArr_JNK3', 2)
Parameter('kr_MKK7BindArr_JNK3', 2)

# Arrestin3-MKK7 bind to uuJNK3, K_D = 1.4 microM, Figure 1.E
Parameter('kf_MKK7_ArrBinduuJNK3', 2)
Parameter('kr_MKK7_ArrBinduuJNK3', 2.4)

# uJNK3 with MKK4
Parameter('kf_MKK4_uJNK3', 2)
Parameter('kr_MKK4_uJNK3', 4)

# uJNK3 with MKK7
# This is when MKK7 is bound
Parameter('kf_MKK7_uJNK3', 2)
Parameter('kr_MKK7_uJNK3', 4)

# Equilibration between unphosphorilated MKK4 and 7, these parameters have units of inverse s
Parameter('keq_uMKK4_to_uMKK7', 10) # 7
Parameter('keq_uMKK7_to_uMKK4', 10) # 4

# MKK4, Arrestin JNK3 activation
Parameter('kcat_pMKK4_ArrJNK3', 1)

# MKK7, Arrestin JNK3 activation
Parameter('kcat_pMKK7_ArrJNK3', 1)

# MKK4, JNK3 activation
Parameter('kcat_pMKK4_JNK3', 1)

# MKK7, JNK3 activation
Parameter('kcat_pMKK7_JNK3', 1)

Parameter('keq_pMKK4_to_pMKK7', 1)
Parameter('keq_pMKK7_to_pMKK4', 1)

Parameter('kf_pJNK3_MKK4complex', 2)
Parameter('kr_pJNK3_MKK4complex', 3)

Parameter('kf_pJNK3_MKK7complex', 2)
Parameter('kr_pJNK3_MKK7complex', 3)

# Initial conditions
Parameter('Arrestin_0', 20)
Parameter('pMKK4_0', 0.2)
Parameter('pMKK7_0', 0.2)
Parameter('uJNK3_0', 2)

Initial(Arrestin(b1=None, b2=None, b3=None), Arrestin_0)
Initial(MKK4(b=None, state='P'), pMKK4_0)
Initial(MKK7(b=None, state='P'), pMKK7_0)
Initial(JNK3(b=None, threo='U', tyro='U'), uJNK3_0)

# Rules

Rule('pMKK4BindArr', Arrestin(b1=None, b2=None, b3=None) + MKK4(b=None, state='P') |
     Arrestin(b1=None, b2=2, b3=None) % MKK4(b=2, state='P'), kf_pMKK4_Arr, kr_pMKK4_Arr)

Rule('pMKK7BindArr', Arrestin(b1=None, b2=None, b3=None) + MKK7(b=None, state='P') |
     Arrestin(b1=None, b2=2, b3=None) % MKK7(b=2, state='P'), kf_pMKK7_Arr, kr_pMKK7_Arr)

Rule('JNK3BindArr', Arrestin(b1=None, b2=None, b3=None) + JNK3(b=None, threo='U', tyro='U') |
     Arrestin(b1=None, b2=None, b3=3) % JNK3(b=3, threo='U', tyro='U'), kf_JNK3_Arr, kr_JNK3_Arr)

Rule('MKK4_ArrBinduuJNK3', Arrestin(b1=None, b2=2, b3=None) % MKK4(b=2, state='P') + JNK3(b=None, threo='U', tyro='U') |
     Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, threo='U', tyro='U'),
     kf_MKK4_ArrBinduuJNK3, kr_MKK4_ArrBinduuJNK3)

# Does the state of JNK3 affect the catalysis?
Rule('MKK4catJNK3Arr', Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, tyro='U') >>
     Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, tyro='P'),kcat_pMKK4_ArrJNK3)

# JNK3 has to get dissociated because otherwise it wouldnt be possible to have more pJNK3 than the value of MKK4 or MKK7
Rule('upJNK3Arr_MKK4_diss', Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3,threo='U', tyro='P') |
     Arrestin(b1=None, b2=2, b3=None) % MKK4(b=2, state='P') + JNK3(b=None, threo='U', tyro='P')
     , kr_upJNK3Arr_MKK4, kf_upJNK3Arr_MKK4)

# Here we assume that the kinetics of puJNK3 disocciation for MKK4 are the same as for MKK7
Rule('puJNK3Arr_MKK4_diss', Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, threo='P', tyro='U') |
     Arrestin(b1=None, b2=2, b3=None) % MKK4(b=2, state='P') + JNK3(b=None, threo='P', tyro='U')
     , kr_puJNK3Arr_MKK7, kf_puJNK3Arr_MKK7)

Rule('ppJNK3Arr_MKK4_diss', Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, threo='P', tyro='P') |
     Arrestin(b1=None, b2=2, b3=None) % MKK4(b=2, state='P') + JNK3(b=None, threo='P', tyro='P')
     , kr_ppJNK3_Arr, kf_ppJNK3_Arr)

Rule('MKK7_ArrBindUUJNK3', Arrestin(b1=None, b2=2, b3=None) % MKK7(b=2, state='P') + JNK3(b=None, threo='U', tyro='U') |
     Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='U', tyro='U')
     , kf_MKK7_ArrBinduuJNK3, kr_MKK7_ArrBinduuJNK3)

Rule('MKK7catJNK3Arr', Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='U') >>
     Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='P'), kcat_pMKK7_ArrJNK3)

# JNK3 has to get dissociated because otherwise it wouldnt be possible to have more pJNK3 than the value of MKK4 or MKK7
Rule('puJNK3Arr_MKK7_diss', Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='P', tyro='U') |
     Arrestin(b1=None, b2=2, b3=None) % MKK7(b=2, state='P') + JNK3(b=None, threo='P', tyro='U')
     , kr_puJNK3Arr_MKK7, kf_puJNK3Arr_MKK7)

# I am not sure about this rule's rate. Does it change the kinetic rate depending if MKK4 or 7 is bound?
# Here we assume that the kinetics of puJNK3 disocciation for MKK4 are the same as for MKK7
Rule('upJNK3Arr_MKK7_diss', Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='U', tyro='P') |
     Arrestin(b1=None, b2=2, b3=None) % MKK7(b=2, state='P') + JNK3(b=None, threo='U', tyro='P')
     , kr_upJNK3Arr_MKK4, kf_upJNK3Arr_MKK4)

Rule('ppJNK3Arr_MKK7_diss', Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='P', tyro='P') |
     Arrestin(b1=None, b2=2, b3=None) % MKK7(b=2, state='P') + JNK3(b=None, threo='P', tyro='P')
     , kr_ppJNK3_Arr, kf_ppJNK3_Arr)

# MKK4/7 release from Arrestin complex
# Does MKK 4/7 bind at a different rate when JNK3 is present? Does it affect if JNK3 is phosphorylated?
Rule('MKK4DissArr_JNK3', Arrestin(b1=None, b2=None, b3=3) % JNK3(b=3) + MKK4(b=None, state='P')|
      Arrestin(b1=None, b2=2, b3=3) % JNK3(b=3) % MKK4(b=2, state='P'), kf_MKK4BindArr_JNK3, kr_MKK4BindArr_JNK3)

Rule('MKK7DissArr_JNK3', Arrestin(b1=None, b2=None, b3=3) % JNK3(b=3) + MKK7(b=None, state='P')|
      Arrestin(b1=None, b2=2, b3=3) % JNK3(b=3) % MKK7(b=2, state='P'), kf_MKK7BindArr_JNK3, kr_MKK7BindArr_JNK3)

# EquilibratePMKK4and7
Rule('EqpMKK4And7', Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, threo='U', tyro='P') + MKK7(b=None, state='P') >>
     Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='U', tyro='P') + MKK4(b=None, state='P'),
     keq_pMKK4_to_pMKK7)

Rule('EqpMKK7And4', Arrestin(b1=None, b2=2, b3=3) % MKK7(b=2, state='P') % JNK3(b=3, threo='P', tyro='U') + MKK4(b=None, state='P') >>
     Arrestin(b1=None, b2=2, b3=3) % MKK4(b=2, state='P') % JNK3(b=3, threo='P', tyro='U') + MKK7(b=None, state='P'),
     keq_pMKK7_to_pMKK4)

#### Interactions that have to be added in order to the model to fit the experimental data

Rule('MKK4BindJNK3', MKK4(b=None, state='P') + JNK3(b=None, tyro='U') |
     MKK4(b=1, state='P') % JNK3(b=1, tyro='U')
     , kf_MKK4_uJNK3, kr_MKK4_uJNK3)

Rule('MKK4catJNK3', MKK4(b=1, state='P') % JNK3(b=1, tyro='U') >>
     MKK4(b=1, state='P') % JNK3(b=1, tyro='P'),kcat_pMKK4_JNK3)

# JNK3 has to get dissociated because otherwise it wouldnt be possible to have more pJNK3 than the value of MKK4 or MKK7
Rule('pJNK3_MKK4complex_diss', MKK4(b=1, state='P') % JNK3(b=1, tyro='P') |
     MKK4(b=None, state='P') + JNK3(b=None, tyro='P'), kr_pJNK3_MKK4complex, kf_pJNK3_MKK4complex)

Rule('MKK7BindJNK3', MKK7(b=None, state='P') + JNK3(b=None, threo='U') |
     MKK7(b=1, state='P') % JNK3(b=1, threo='U'), kf_MKK7_uJNK3, kr_MKK7_uJNK3)

Rule('MKK7catJNK3', MKK7(b=1, state='P') % JNK3(b=1, threo='U') >>
     MKK7(b=1, state='P') % JNK3(b=1, threo='P'), kcat_pMKK7_JNK3)

# JNK3 has to get dissociated because otherwise it wouldnt be possible to have more pJNK3 than the value of MKK4 or MKK7
Rule('pJNK3_MKK7complex_diss', MKK7(b=1, state='P') % JNK3(b=1, threo='P') |
     MKK7(b=None, state='P') + JNK3(b=None, threo='P') , kr_pJNK3_MKK7complex, kf_pJNK3_MKK7complex)

# Unbound JNK3
Observable('pTyr_jnk3', JNK3(b=None, tyro='P'))
Observable('pThr_jnk3', JNK3(b=None, threo='P'))
