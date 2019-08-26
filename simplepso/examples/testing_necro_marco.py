
from pysb import *
import pandas as pd
from pysb.core import *
from pysb.bng import *
from pysb.integrate import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import *
import scipy.interpolate
import numpy as np
from pysb.simulator.bng import BngSimulator

Model ()

#Declaration of monomers

""" Delcares IkBa, IkBb, IkBe, and IkBd monomers. All four of these monomers have binding sites to IKK2
    (for IkBa, IkBb, IkBe), and IKK1 (for IkBd). The monomers import and export between the Cytoplasm
    and the Nucleus.  Declares the IkBa, IkBb, IkBe, and IkBd mRNA monomers. The monomers go through constitutive RNA synthesis
    without the presence of NFkB in the nucleus,  and inducible RNA synthesis in the presence of NFkB in the nucleus. """

Monomer('IkBa', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBa_mRNA')
Monomer('IkBb', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBb_mRNA')
Monomer('IkBe', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBe_mRNA')
Monomer('IkBd', ['nfkb', 'S'], {'S': ['C', 'N']})
Monomer('IkBd_mRNA')
Monomer('gene', ['type','state'], {'type': ['a', 'b', 'e', 'd'],'state':['off', 'on']})
Monomer('NFkB', ['ikb', 'S'], {'S': ['C', 'N']})
Monomer('A20')
Monomer('A20t')
Monomer('TNF', ['c', 'tnfr'])
Monomer('TNFRM')
Monomer('TNFR', ['tnf'])
Monomer('C1', ['tnf', 'state'], {'state': ['a', 'i']})
Monomer('TTR')
Monomer('IKKK', ['state'], {'state': ['a', 'i']})
Monomer('IKK', ['state'], {'state': ['a', 'i', 'ai']})

Parameter('IkBa_0', 0.0025) #IkBa
Initial(IkBa(nfkb = None, S = 'C'), IkBa_0)
Parameter('IkBan_0', 0.0013) #IkBa
Initial(IkBa(nfkb = None, S = 'N'), IkBan_0)
Parameter('IkBa_NFkB_0', 0.0621) #IkBa
Initial(IkBa(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBa_NFkB_0)
Parameter('IkBa_NFkBn_0', 0.0208) #IkBa
Initial(IkBa(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBa_NFkBn_0)
Parameter('IkBat_0', 0.0020) #IkBa
Initial(IkBa_mRNA(), IkBat_0)
Parameter('IkBb_0', 0.0044) #IkBb
Initial(IkBb(nfkb = None, S = 'C'), IkBb_0)
Parameter('IkBbn_0', 0.0002) #IkBa
Initial(IkBb(nfkb = None, S = 'N'), IkBbn_0)
Parameter('IkBb_NFkB_0', 0.0237) #IkBa
Initial(IkBb(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBb_NFkB_0)
Parameter('IkBb_NFkBn_0', 0.0016) #IkBb
Initial(IkBb(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBb_NFkBn_0)
Parameter('IkBbt_0', 0.0033) #IkBa
Initial(IkBb_mRNA(), IkBbt_0)
Parameter('IkBe_0', 0.0003) #IkBa
Initial(IkBe(nfkb = None, S = 'C'), IkBe_0)
Parameter('IkBen_0', 0.0001) #IkBa
Initial(IkBe(nfkb = None, S = 'N'), IkBen_0)
Parameter('IkBeNFkB_0', 0.0045) #IkBa
Initial(IkBe(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBeNFkB_0)
Parameter('IkBeNFkBn_0', 0.0015) #IkBa
Initial(IkBe(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBeNFkBn_0)
Parameter('IkBet_0', 0.0003) #IkBa
Initial(IkBe_mRNA(), IkBet_0)
Parameter('IkBd_0', 0.0003) #IkBa
Initial(IkBd(nfkb = None, S = 'C'), IkBd_0)
Parameter('IkBdn_0', 0.0003) #IkBa
Initial(IkBd(nfkb = None, S = 'N'), IkBdn_0)
Parameter('IkBdNFkB_0', 0.0058) #IkBa
Initial(IkBd(nfkb = 1, S = 'C')%NFkB(ikb = 1, S = 'C'), IkBdNFkB_0)
Parameter('IkBdNFkBn_0', 0.0039) #IkBa
Initial(IkBd(nfkb = 1, S = 'N')%NFkB(ikb = 1, S = 'N'), IkBdNFkBn_0)
Parameter('IkBdt_0', 0.0001) #IkBa
Initial(IkBd_mRNA(), IkBdt_0)
Parameter('NFkB_0', 0.000) #Nuclear Factor-kappaB
Initial(NFkB(ikb=None, S='C'), NFkB_0)
Parameter('NFkBn_0', 0.0012) #Nuclear Factor-kappaB in nucleus
Initial(NFkB(ikb=None, S='N'), NFkBn_0)
Parameter('IKKK_off_0', 0.1000) #Inhibitor Kinase Kinase Kinase
Initial(IKKK(state = 'i'), IKKK_off_0)
Parameter('IKKK_0', 0.000) #Inhibitor Kinase Kinase Kinase
Initial(IKKK(state = 'a'), IKKK_0)
Parameter('IKK_off_0', 0.0985) #TRADD-TRAF-RIP
Initial(IKK(state = 'i'), IKK_off_0)
Parameter('IKK_0', 0.0002) #TRADD-TRAF-RIP
Initial(IKK(state = 'a'), IKK_0)
Parameter('IKKi_0', 0.0013) #TRADD-TRAF-RIP
Initial(IKK(state = 'ai'), IKKi_0)
Parameter('TNF_0', .2) #TNF
Initial(TNF(c = None, tnfr = None), TNF_0)
Parameter('TNFRM_0', 0.0) #TRADD-TRAF-RIP
Initial(TNFRM(), TNFRM_0)
Parameter('TNFR_0', 0.0) #TNFR Tumor Necrosis Factor Receptor 1
Initial(TNFR(tnf = None), TNFR_0)
Parameter('TNFRtnf_0', 0.0) #TNF binding TNFR
Initial(TNF(c = None, tnfr = 1)%TNFR(tnf = 1), TNFRtnf_0)
Parameter('C1_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = None, state = 'a'), C1_0)
Parameter('C1_off_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = None, state = 'i'), C1_off_0)
Parameter('C1tnf_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = 1, state = 'a')%TNF(c =  1, tnfr = None), C1tnf_0)
Parameter('C1tnf_off_0', 0.0) #TRADD-TRAF-RIP
Initial(C1(tnf = 1, state = 'i')%TNF(c =  1, tnfr = None), C1tnf_off_0)
Parameter('TTR_0', 8.3e-4) #TRADD-TRAF-RIP
Initial(TTR(), TTR_0)
Parameter('a20_0', 0.0049) #TRADD-TRAF-RIP
Initial(A20(), a20_0)
Parameter('a20t_0', 0.0001) #TRADD-TRAF-RIP
Initial(A20t(), a20t_0)
Parameter('genea_0', 0.012)
Initial(gene(type = 'a', state = 'off'), genea_0)
Parameter('geneb_0', 0.012)
Initial(gene(type = 'b', state = 'off'), geneb_0)
Parameter('genee_0', 0.012)
Initial(gene(type = 'e', state = 'off'), genee_0)
Parameter('gened_0', 0.012)
Initial(gene(type = 'd', state = 'off'), gened_0)

# def  ikb_mrna_to_ikb():
Parameter('psynth_a', 7e-5)
Parameter('psynth_b', 1e-5)
Parameter('psynth_e', 1e-6)
Parameter('psynth_d', 1e-7)
Rule('a_synth', None >> IkBa_mRNA(), psynth_a)
Rule('b_synth', None >> IkBb_mRNA(), psynth_b)
Rule('e_synth', None >> IkBe_mRNA(), psynth_e)
Rule('d_synth', None >> IkBd_mRNA(), psynth_d)

Parameter('a_delay', 0.1)
Parameter('be_delay', 0.027)
Parameter('d_delay', 0.011)
Rule('genea_off_on', gene(type = 'a', state = 'off') >> gene(type = 'a',state = 'on'), a_delay)
Rule('geneb_off_on', gene(type = 'b', state = 'off') >> gene(type = 'b',state = 'on'), be_delay)
Rule('genee_off_on', gene(type = 'e', state = 'off') >> gene(type = 'e',state = 'on'), be_delay)
Rule('gened_off_on', gene(type = 'd', state = 'off') >> gene(type = 'd',state = 'on'), d_delay)

Parameter('hill', 3)
Parameter('a', 8)
Parameter('b', 0.02)
Parameter('e', 0.3)
Parameter('d', 0.025)
Observable('NFkBn_free', NFkB(ikb=None, S='N'))
Expression('a_NFkBn', a*(NFkBn_free)**(hill))
Expression('b_NFkBn', b*(NFkBn_free)**(hill))
Expression('e_NFkBn', e*(NFkBn_free)**(hill))
Expression('d_NFkBn', d*(NFkBn_free)**(hill))
Rule('an_mRNA', gene(type = 'a',state = 'on') >> IkBa_mRNA() + gene(type = 'a',state = 'on'), a_NFkBn)
Rule('bn_mRNA', gene(type = 'b',state = 'on') >> IkBb_mRNA() + gene(type = 'b',state = 'on'), a_NFkBn)
Rule('en_mRNA', gene(type = 'e',state = 'on') >> IkBe_mRNA() + gene(type = 'e',state = 'on'), a_NFkBn)
Rule('dn_mRNA', gene(type = 'd',state = 'on') >> IkBd_mRNA() + gene(type = 'd',state = 'on'), a_NFkBn)

Parameter('mRNA_a', 0.035)
Parameter('mRNA_b', 3e-3)
Parameter('mRNA_e', 4e-3)
Parameter('mRNA_d', 2e-3)
Rule('a_mRNA', IkBa_mRNA() >> None, mRNA_a)
Rule('b_mRNA', IkBb_mRNA() >> None, mRNA_b)
Rule('e_mRNA', IkBe_mRNA() >> None, mRNA_e)
Rule('d_mRNA', IkBd_mRNA() >> None, mRNA_d)


Parameter('syntha', 0.2448)
Parameter('synthb', 0.2448)
Parameter('synthe', 0.2448)
Parameter('synthd', 0.2448)
Rule('a_psynth', IkBa_mRNA() >> IkBa(nfkb=None, S='C') + IkBa_mRNA(), syntha)
Rule('b_psynth', IkBb_mRNA() >> IkBb(nfkb=None, S='C') + IkBb_mRNA(), synthb)
Rule('e_psynth', IkBe_mRNA() >> IkBe(nfkb=None, S='C') + IkBe_mRNA(), synthe)
Rule('d_psynth', IkBd_mRNA() >> IkBd(nfkb=None, S='C') + IkBd_mRNA(), synthd)


Parameter('IkB_IKKf', 30)
Parameter('IkB_IKKr', 6e-5)
Rule('an_adc', IkBa(nfkb=None, S='C') + NFkB(ikb=None, S='C') | IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('bn_adc', IkBb(nfkb=None, S='C') + NFkB(ikb=None, S='C') | IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('en_adc', IkBe(nfkb=None, S='C') + NFkB(ikb=None, S='C') | IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)
Rule('dn_adc', IkBd(nfkb=None, S='C') + NFkB(ikb=None, S='C') | IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C'), IkB_IKKf, IkB_IKKr)

Rule('an_adn', IkBa(nfkb=None, S='N') + NFkB(ikb=None, S='N') | IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('bn_adn', IkBb(nfkb=None, S='N') + NFkB(ikb=None, S='N') | IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('en_adn', IkBe(nfkb=None, S='N') + NFkB(ikb=None, S='N') | IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)
Rule('dn_adn', IkBd(nfkb=None, S='N') + NFkB(ikb=None, S='N') | IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'), IkB_IKKf, IkB_IKKr)


Parameter('af', 0.09)
Parameter('bf', 0.009)
Parameter('ef', 0.045)
Parameter('df', 0.045)
Parameter('ancf', 0.012)
Parameter('bncf', 0.012)
Parameter('encf', 0.012)
Parameter('dncf', 0.012)
Rule('a_nc', IkBa(nfkb=None, S='C') | IkBa(nfkb=None, S='N'), af, ancf)
Rule('b_nc', IkBb(nfkb=None, S='C') | IkBb(nfkb=None, S='N'), bf, bncf)
Rule('e_nc', IkBe(nfkb=None, S='C') | IkBe(nfkb=None, S='N'), ef, encf)
Rule('d_nc', IkBd(nfkb=None, S='C') | IkBd(nfkb=None, S='N'), ef, dncf)

Parameter('anf', 0.276)
Parameter('bnf', 0.0276)
Parameter('enf', 0.138)
Parameter('dnf', 0.276)
Parameter('anr', 0.828)
Parameter('bnr', 0.414)
Parameter('enr', 0.414)
Parameter('dnr', 0.414)
Rule('an_nc', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') | IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N'), anf, anr)
Rule('bn_nc', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') | IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N'), bnf, bnr)
Rule('en_nc', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') | IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N'), enf, enr)
Rule('dn_nc', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') | IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N'), dnf, dnr)

Parameter('nf', 5.4)
Parameter('nr', 0.0048)
Rule('n_nc', NFkB(ikb=None, S='C') | NFkB(ikb=None, S='N'), nf, nr)

Parameter('a_dc', 0.12)
Parameter('b_dc', 0.18)
Parameter('e_dc', 0.18)
Parameter('d_dc', 0.0014)
Parameter('a_dn', 0.12)
Parameter('b_dn', 0.18)
Parameter('e_dn', 0.18)
Parameter('d_dn', 0.0014)

Rule('ad_c', IkBa(nfkb=None, S='C') >> None, a_dc)
Rule('bd_c', IkBb(nfkb=None, S='C') >> None, b_dc)
Rule('ed_c', IkBe(nfkb=None, S='C') >> None, e_dc)
Rule('dd_c', IkBd(nfkb=None, S='C') >> None, d_dc)

Rule('ad_n', IkBa(nfkb=None, S='N') >> None, a_dn)
Rule('bd_n', IkBb(nfkb=None, S='N') >> None, b_dn)
Rule('ed_n', IkBe(nfkb=None, S='N') >> None, e_dn)
Rule('dd_n', IkBd(nfkb=None, S='N') >> None, d_dn)

Parameter('c_bn', 0.00006)
Parameter('n_bn', 0.00006)

Rule('an_c', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('bn_c', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('en_c', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)
Rule('dn_c', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), c_bn)

Rule('an_n', IkBa(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('bn_n', IkBb(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('en_n', IkBe(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)
Rule('dn_n', IkBd(nfkb=1, S='N')%NFkB(ikb=1, S='N') >> NFkB(ikb=None, S='N'), n_bn)



Parameter('a_f_deg', 0.36/.1)
Parameter('b_f_deg', 0.12/.1)
Parameter('e_f_deg', 0.18/.1)
Parameter('d_f_deg', 0.18*.01)

# Parameter('and_c_n', 0.36)
# Parameter('bnd_c_n', 0.12)
# Parameter('end_c_n', 0.18)
Parameter('dnd_c_n', 0.18*.01)

Observable('IKKa_obs', IKK(state = 'a'))

Expression('IKK_ikba_flux', model.observables['IKKa_obs']*a_f_deg)
Expression('IKK_ikbb_flux', model.observables['IKKa_obs']*b_f_deg)
Expression('IKK_ikbe_flux', model.observables['IKKa_obs']*e_f_deg)

Rule('a_c_deg', IkBa(nfkb=None, S='C') >> None, IKK_ikba_flux)
Rule('b_c_deg', IkBb(nfkb=None, S='C') >> None, IKK_ikbb_flux)
Rule('e_c_deg', IkBe(nfkb=None, S='C') >> None, IKK_ikbe_flux)
Rule('d_c_deg', IkBd(nfkb=None, S='C') >> None, d_f_deg)

Rule('an_c_n', IkBa(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikba_flux)
Rule('bn_c_n', IkBb(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikbb_flux)
Rule('en_c_n', IkBe(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), IKK_ikbe_flux)
Rule('dn_c_n', IkBd(nfkb=1, S='C')%NFkB(ikb=1, S='C') >> NFkB(ikb=None, S='C'), dnd_c_n)

Parameter('A20_mRNA', 2e-6)
Parameter('A20n', 0.4)
Parameter('A20_mRNA_c_deg', 0.035)
# Parameter('a1d_c_deg', 0.36)
Parameter('A20_synth', 0.25)
Parameter('A20_deg', 0.0029)

Observable('obs_A20t', A20t())

Expression('A20t_NFkBn', A20n*(NFkBn_free)**(hill))
Expression('A20_synthesis', A20_synth*obs_A20t)

Rule('A20t_synth', None >> A20t(), A20_mRNA)
Rule('A20t_mediated_nfkbn', None >> A20t(), A20t_NFkBn)
Rule('A20t_deg', A20t() >> None, A20_mRNA_c_deg)

Rule('synth_A20', None >> A20(), A20_synthesis)
Rule('deg_A20', A20() >> None, A20_deg)

Parameter('synth_tnfrm', 2e-7)
Parameter('deg_tnfrm', 0.0058)
Parameter('tnfr_f_tnfrm', 1e-5*3) # times 2
Parameter('tnfr_r_tnfrm', 0.1*3)
Parameter('deg_TNFR', 0.023)
Parameter('TNFR_TTR_f_C1', 100.0)
Parameter('TNFR_TTR_r_C1', 0.1)
Parameter('C1_f_a', 30.0)
Parameter('C1_r_a', 2.0)
Parameter('C1_f_A20', 1000.0)
Parameter('C1_f_TNFR_TTR', 0.1)
Parameter('C1_i_deg', 0.023)
Parameter('C1_a_deg', 0.023)
Rule('tnfrm_synth', None >> TNFRM(), synth_tnfrm)
Rule('tnfrm_deg', TNFRM() >> None, deg_tnfrm)
Rule('tnfrm3_to_TNFR', TNFRM()  >> TNFR(tnf = None), tnfr_f_tnfrm)
Rule('TNFR_to_tnfrm3', TNFR(tnf = None) >> TNFRM() , tnfr_r_tnfrm)
Rule('TNFR_deg', TNFR(tnf = None) >> None, deg_TNFR)
Rule('TNFR_TTR_C1i', TNFR(tnf = None) + TTR() | C1(tnf = None, state = 'i'), TNFR_TTR_f_C1, TNFR_TTR_r_C1)
Rule('a_C1_i', C1(tnf = None, state = 'i') | C1(tnf = None, state = 'a'), C1_f_a, C1_r_a)
Rule('C1_i_A20', C1(tnf = None, state = 'a') + A20() >> C1(tnf = None, state = 'i') + A20(), C1_f_A20)
Rule('C1_TNFR_TTR', C1(tnf = None, state = 'a') >> TNFR(tnf = None) + TTR(), C1_f_TNFR_TTR)
Rule('C1i_deg', C1(tnf = None, state = 'i') >> None, C1_i_deg)
Rule('C1a_deg', C1(tnf = None, state = 'a') >> None, C1_a_deg)

Parameter('tnf_deg', 0.0154)
Parameter('tnf_tnfrm_f_TNFRtnf', 1100.0*3) # times 2
Parameter('tnf_TNFR_f_TNFRtnf', 1100.0)
Parameter('tnf_TNFR_r_TNFRtnf', 0.021)
Parameter('deg_TNFRtnf', 0.023)
Parameter('TNFRtnf_TTR_f_C1itnf', 100)
Parameter('TNFRtnf_TTR_r_C1itnf', 0.1)
Parameter('C1itnf_f_C1atnf', 30.0)
Parameter('C1itnf_r_C1atnf', 2.0)
Rule('deg_tnf', TNF(c = None, tnfr = None) >> None, tnf_deg)
Rule('tnf_tnfrm_TNFRtnf', TNF(c = None, tnfr = None) + TNFRM()  >> TNFR(tnf = 1)%TNF(c=None, tnfr = 1), tnf_tnfrm_f_TNFRtnf)
Rule('tnf_TNFR_TNFRtnf', TNFR(tnf = None) + TNF(c = None, tnfr = None)  | TNFR(tnf = 1)%TNF(c=None, tnfr = 1), tnf_TNFR_f_TNFRtnf, tnf_TNFR_r_TNFRtnf)
Rule('TNFRtnf_deg', TNFR(tnf = 1)%TNF(c =None, tnfr = 1) >> None, deg_TNFRtnf)
Rule ('TNFRtnf_TTR_C1itnf', TNFR(tnf = 1) % TNF(c = None, tnfr = 1) + TTR() | C1(tnf = 1, state = 'i') % TNF(c = 1, tnfr =None), TNFRtnf_TTR_f_C1itnf, TNFRtnf_TTR_r_C1itnf)
Rule('C1itnf_C1atnf', C1(tnf = 1, state = 'i') % TNF(c = 1, tnfr = None) | C1(tnf = 1, state = 'a') % TNF(c = 1, tnfr = None), C1itnf_f_C1atnf, C1itnf_r_C1atnf)

Parameter('C1atnf_f_C1itnf', 1000.0)
Parameter('C1tnf_f_TNFRtnf_TTR', 0.1)
Parameter('deg_C1itnf', 0.023)
Parameter('deg_C1atnf', 0.023)
Rule('C1atnf_C1itnf',C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) + A20() >> C1(tnf = 1, state = 'i')%TNF(c = 1, tnfr = None) + A20(), C1atnf_f_C1itnf)
Rule('C1tnf_TNFRtnf_TTR', C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) >> TNFR(tnf = 1)%TNF(c = None, tnfr = 1) + TTR(), C1tnf_f_TNFRtnf_TTR)
Rule('C1itnf_deg', C1(tnf = 1, state = 'i')%TNF(c = 1, tnfr = None) >> None, deg_C1itnf)
Rule('C1atnf_deg', C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) >> None, deg_C1atnf)

Parameter('C1i_f_tnf', 0.021)
Parameter('C1i_r_tnf', 1100.0)
Parameter('C1a_f_tnf', 0.021)
Parameter('C1a_r_tnf', 1100.0)
Rule('C1itnf', C1(tnf = 1, state = 'i')%TNF(c = 1, tnfr = None) | C1(tnf = None, state = 'i') + TNF(c = None, tnfr = None), C1i_f_tnf, C1i_r_tnf)
Rule('C1atnf', C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) | C1(tnf = None, state = 'a') + TNF(c = None, tnfr = None), C1a_f_tnf, C1a_r_tnf)

Parameter('IKKKi_f_IKKK', 5e-7)
Parameter('IKKKa_r_IKKK', 0.25)
Parameter('IKKKi_f_IKKKaC1', 500.0)
Parameter('IKKKi_f_IKKKaC1tnf', 500.0)
Rule('IKKKi_IKKKa', IKKK(state = 'i') | IKKK(state = 'a'), IKKKi_f_IKKK, IKKKa_r_IKKK)
Rule('IKKKi_IKKKaC1', IKKK(state = 'i') + C1(tnf = None, state = 'a') >> IKKK(state = 'a') + C1(tnf = None, state = 'a'), IKKKi_f_IKKKaC1)
Rule('IKKKi_IKKKaC1tnf', IKKK(state = 'i') + C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None) >> IKKK(state = 'a') + C1(tnf = 1, state = 'a')%TNF(c = 1, tnfr = None), IKKKi_f_IKKKaC1tnf)

#IKK Activity Reactions
Parameter('IKKi_f_IKK', 5e-5)
Parameter('IKKa_r_IKK', 0.02)
Parameter('IKKi_f_IKKaIKKK', 520.0)
Parameter('IKKa_f_IKKai', 0.15)
Parameter('IKKai_f_IKKi', 0.02)
Rule('IKKi_IKKa', IKK(state = 'i') | IKK(state = 'a'), IKKi_f_IKK, IKKa_r_IKK)
Rule('IKKi_IKKaIKKK', IKK(state = 'i') + IKKK(state = 'a') >> IKK(state = 'a') + IKKK(state = 'a'), IKKi_f_IKKaIKKK)
Rule('IKKi_IKKai', IKK(state = 'a') >> IKK(state = 'ai'), IKKa_f_IKKai)
Rule('IKKai_IKKi', IKK(state = 'ai') >> IKK(state = 'i'), IKKai_f_IKKi)

Observable('NFkBn_obs', NFkB(ikb=None, S='N'))
Observable('TNFRM_obs', TNFRM())
Observable('TNFR_obs', TNFR(tnf = None))
Observable('IkBa_mRNA_obs', IkBa_mRNA())
Observable('IkBb_mRNA_obs', IkBb_mRNA())
Observable('IkBe_mRNA_obs', IkBe_mRNA())
Observable('IkBd_mRNA_obs', IkBd_mRNA())
Observable('TNF_obs', TNF(c = None, tnfr = None))
Observable('IkBa_obs', IkBa(nfkb = None, S='C'))
Observable('IkBb_obs', IkBb(nfkb = None, S='C'))
Observable('IkBe_obs', IkBe(nfkb = None, S='C'))
Observable('IkBd_obs', IkBd(nfkb = None, S='C'))

generate_equations(model)
# data = np.loadtxt('../output_population_nfkb')
# print(data.shape)
print(len(model.parameters_rules()))
print(len(model.parameters_expressions()))
print(len(model.initial_conditions))
print(model.parameters)
print(len(model.odes))
print(len(model.species))
quit()

tspan = np.linspace(0, 720, 721)
sim = BngSimulator(model, tspan = tspan)
simulation_result = sim.run(method='ode')

plt.figure(figsize = (15,10))
# plt.figure()
plt.subplot(231)
plt.plot(tspan/60, simulation_result.observables['TNF_obs'], label = 'TNF')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(232)
plt.plot(tspan/60, simulation_result.observables['TNFR_obs'], label = 'TNFR')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(233)
plt.plot(tspan/60, simulation_result.observables['IKKa_obs'],label = 'IKKa')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

plt.subplot(234)
plt.plot(tspan/60, simulation_result.observables['NFkBn_obs'], label = 'NFkBn')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(235)
plt.plot(tspan/60, simulation_result.observables['obs_A20t'], label = 'A20')

plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)

# plt.figure()
plt.subplot(236)
plt.plot(tspan/60, simulation_result.observables['IkBa_obs'],label = 'IkBa')
plt.xlabel("Time (in hr)", fontsize=10)
plt.ylabel("Concentration uM", fontsize=10)
# plt.ylim(ymin = -10, ymax =100)
plt.legend(loc=0)
plt.savefig('bng_nfkb')
plt.tight_layout()
plt.show()