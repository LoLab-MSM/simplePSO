# MW Irvin -- Lopez Lab -- 2018-10-01
"""
========================================
Simple PySB Model of Apoptosis Signaling
========================================
PySB Model of the following Apoptosis signaling reactions:
Bid Independent Apoptosis
1. Receptor ligation:                    TRAIL + Receptor  <-> TRAIL:Receptor     kf0, kr0
2. Death inducing signaling complex:       TRAIL:Receptor  --> DISC               kc0
3. Initiator Caspase Activation: DISC->         iCaspases  --> *iCaspases         kf1, kr1, kc1
4. Effector Caspase Activation:  *iCaspases->   eCaspases  --> *eCaspases         kf2, kr2, kc2
5. Feedback Capsapse Activation: *eCaspases->   iCaspases  --> *iCaspases         kf3, kr3, kc3
Bid Dependent Apoptosis
6. Bid Activation:               *iCaspases->          Bid --> *Bid               kf4, kr4, kc4
7. MOMP Dependent Signaling:           *Bid-> MOMP signals --> *MOMP signals      kf5, kr5, kc5
8. Effector Caspase Activation: *MOMP signals -> eCaspases --> *eCaspases         kf6, kr6, kc6
9. PARP Cleavage:                  *eCaspases->  PARP      --> cPARP              kf7, kr7, kc7
10. Initiator Caspase Degradation:              iCaspsases --> None               kc8
"""
from pysb import *
from pysb.macros import catalyze_state, degrade
Model()

# Model reaction rates have units of copies per cell per min

Parameter('L_0',       3000.0)   # baseline level of TRAIL in most experiments (50 ng/ml SuperKiller TRAIL)
Parameter('R_0',        200.0)   # TRAIL receptor (for experiments not involving siRNA)
Parameter('DISC_0',       0.0)   # Death inducing signaling complex
Parameter('IC_0',       2.0e4)   # Initiator Caspases
Parameter('EC_0',       1.0e4)   # Effector Caspases
Parameter('Bid_0',      4.0e4)   # Bid
Parameter('MOMP_sig_0', 1.0e5)   # MOMP dependent pro-apoptotic Effector Caspase activators (e.g. cytochrome C)
Parameter('PARP_0',     1.0e6)   # PARP (Caspase-3 substrate)

Monomer('L',          ['b'])
Monomer('R',          ['b'])
Monomer('DISC',       ['b'])
Monomer('IC',         ['b', 'state'], {'state': ['inactive', 'active']})
Monomer('EC',         ['b', 'state'], {'state': ['inactive', 'active']})
Monomer('Bid',        ['b', 'state'], {'state': ['unmod', 'cleaved']})
Monomer('MOMP_sig',   ['b', 'state'], {'state': ['inactive', 'active']})
Monomer('PARP',       ['b', 'state'], {'state': ['unmod', 'cleaved']})

Initial(L(b=None), L_0)
Initial(R(b=None), R_0)
Initial(DISC(b=None), DISC_0)
Initial(IC(b=None, state='inactive'), IC_0)
Initial(EC(b=None, state='inactive'), EC_0)
Initial(Bid(b=None, state='unmod'), Bid_0)
Initial(MOMP_sig(b=None, state='inactive'), MOMP_sig_0)
Initial(PARP(b=None, state='unmod'), PARP_0)

# 1. Receptor ligation: TRAIL + Receptor  <-> TRAIL:Receptor     kf0, kr0
Parameter('kf0', 4.0e-07)
Parameter('kr0', 1.0e-03)
Rule('Receptor_ligation', L(b=None) + R(b=None) | L(b=1) % R(b=1), kf0, kr0)

# 2. Death inducing signaling complex: TRAIL:Receptor  --> DISC  kc0
Parameter('kc0', 1.0e-05)
Rule('DISC_formation', L(b=1) % R(b=1) >> DISC(b=None), kc0)

# 3. Initiator Caspase Activation: DISC->  iCaspases  --> *iCaspases   kf1, kr1, kc1
Parameter('kf1', 1.0e-08)
Parameter('kr1', 1.0e-03)
Parameter('kc1', 1.0e-00)
catalyze_state(DISC(), 'b', IC(), 'b', 'state', 'inactive', 'active', [kf1, kr1, kc1])

# 4. Effector Caspase Activation:  *iCaspases->   eCaspases  --> *eCaspases kf2, kr2, kc2
Parameter('kf2', 1.0e-10)
Parameter('kr2', 1.0e-03)
Parameter('kc2', 1.0e-02)
catalyze_state(IC(state='active'), 'b', EC(), 'b', 'state', 'inactive', 'active', [kf2, kr2, kc2])

# 5. Feedback Capsase Activation: *eCaspases->   iCaspases  --> *iCaspases kf3, kr3, kc3
Parameter('kf3', 3.0e-08)
Parameter('kr3', 1.0e-03)
Parameter('kc3', 1.0e-02)
catalyze_state(EC(state='active'), 'b', IC(), 'b', 'state', 'inactive', 'active', [kf3, kr3, kc3])

# 6. Bid Activation: *iCaspases->  Bid --> *Bid  kf4, kr4, kc4
Parameter('kf4', 1.0e-06)
Parameter('kr4', 1.0e-03)
Parameter('kc4', 1.0e-00)
catalyze_state(IC(state='active'), 'b', Bid(), 'b', 'state', 'unmod', 'cleaved', [kf4, kr4, kc4])

# 7. MOMP Dependent Signaling: *Bid-> MOMP signals --> *MOMP signals      kf5, kr5, kc5
Parameter('kf5', 1.0e-06)
Parameter('kr5', 1.0e-03)
Parameter('kc5', 1.0e-05)
catalyze_state(Bid(state='cleaved'), 'b', MOMP_sig(), 'b', 'state', 'inactive', 'active', [kf5, kr5, kc5])

# 8. Effector Caspase Activation: *MOMP signals -> eCaspases --> *eCaspases  kf6, kr6, kc6
Parameter('kf6', 1.0e-06)
Parameter('kr6', 1.0e-03)
Parameter('kc6', 1.0e-00)
catalyze_state(MOMP_sig(state='active'), 'b', EC(), 'b', 'state', 'inactive', 'active', [kf6, kr6, kc6])

# 9. PARP Cleavage:                *eCaspases->  PARP      --> cPARP              kf7, kr7, kc7
Parameter('kf7', 1.0e-06)
Parameter('kr7', 1.0e-02)
Parameter('kc7', 1.0e-00)
catalyze_state(EC(state='active'), 'b', PARP(), 'b', 'state', 'unmod', 'cleaved', [kf7, kr7, kc7])

# 10. Initiator Caspase Degradation
Parameter('kc8', 1.0e-06)
degrade(IC(state='active'), kc8)

Observable('TRAIL_receptor_obs', L(b=1) % R(b=1) + DISC())
Observable('DISC_obs', DISC())
Observable('C8_DISC_recruitment_obs', DISC()%IC())
Observable('C8_active_obs', IC(state='active'))
Observable('C8_inactive_obs', IC(state='inactive'))
Observable('C3_active_obs', EC(state='active'))
Observable('C3_inactive_obs', EC(state='inactive'))
Observable('BID_obs', Bid(state='unmod'))
Observable('tBID_obs',   Bid(state='cleaved'))
Observable('cPARP_obs', PARP(b=None, state='cleaved'))
Observable('PARP_obs',  PARP(b=None, state='unmod'))