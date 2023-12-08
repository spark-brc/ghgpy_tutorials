
class DCparms(object):
    
    def __init__(self):
        self.cvfr_cho_to_ch4 = 0.5 # conversion factor of carbohydrate decomposition to CO2 and CH4 (dimensionless) 
                                    # C6H12O6_to_CH4 (alpha1)
        self.cvcf_oc_to_co2 = 0.5 # conversion coefficient on a mole weight basis of organic carbon to CO2 
                                    #calculated as 0.5 (alpha2)
        self.frToCH4 = 0.15 # the fraction of heterotrophic respiration converted to CH4 under anaerobic conditions
                            # beta1
        self.frToExudate = 0.45 # fraction of root carbon production contributing to rhizodeposition (fracToExudates)
                            # range 0.35-0.6 beta2 FREXUD 0.0 - 1.0
        self.hr = 0.1 # Rh is heterotrophic respiration from decomposition of organic matter 
                        #(above- and below-ground structural and metabolic litter and above- and below-ground SOC pools) 
                        # (g CO2-C m−2 d−1)
                        # NOTE: is this daily simulated input or parameter?
        self.aeh = 0.23 # range: 
        self.deh = 0.16 # range:
        self.beh_flood = -250 # Lower limit value for Eh during flooding course (mv)
        self.beh_drain = 300
        self.eh_rain_irr = -20
        self.zero_root_frac = 0.7 # (0-1)
        self.ch4rootlim = 1.0 
        self.co2_to_ch4 = 0.5
        self.mxch4f = 0.55 # MaXimum Fraction of CH4 production emitted by plants. (MXCH4F)
        self.frac_to_exd = 0.45 # range(0 - 1)
        self.tmxbio = 1260 # (rice)
        self.q10 = 3.0
        self.frCH4emit_b = 0.7 # the fraction of CH4 emitted via bubbles
        self.mo = 0.002 # Mo was set to 0.0015 gC m−2d−1 (Huang et al. 2004) but is set to 0.002 in DayCent
        self.mrtblm = 1.0 # Root biomass that when exceeded starts to reduce methane bubble formation (g biomass m-2)


class DNDCparms(object):
    """_summary_

    :param object: _description_
    :type object: _type_
    """
    def __init__(
            self, 
            a=1,
            b=1,
            rconst=8.314,
            fconst=96485,
            fr_rsi=0.5,
            fr_root=0.1
            ):
        self.a = a 
        self.b = b  
        self.rconst = rconst  
        self.fconst = fconst  
        self.fr_rsi = fr_rsi
        self.fr_root = fr_root

class MERESparms(object):
    """_summary_

    :param c_aex: the critical concentration of the oxidized alternative electron acceptor pool (mol Ceq m-3), defaults to 24.0
    :type c_aex: float, optional, mol Ceq m-3
    :param eta: representing the sensitivity of methanogenesis to the concentration of O2
                A value of 400 m3 mol-1 was used for η (Arah & Stephen, 1998)., defaults to 400
    :type eta: int, optional, m3 mol-1
    :param k1: k1, k2 and k3 are Michaelis-Menten constants for a dual-substrate reaction., defaults to 0.33
    :type k1: float, optional, mol m-3
    :param k2: k1, k2 and k3 are Michaelis-Menten constants for a dual-substrate reaction., defaults to 0.44
    :type k2: float, optional, mol m-3
    :param k3: k1, k2 and k3 are Michaelis-Menten constants for a dual-substrate reaction., defaults to 0.22
    :type k3: float, optional, mol m-3
    :param sc_root: the specific conductivity of the root system., defaults to 3.0e-04
    :type sc_root: float, optional, m air (m root)-1
    :param o2conc_atm: the concentration of the O2 in the atmosphere, defaults to 7.76
    :type o2conc_atm: float, optional, mol m-3
    :param ch4conc_atm: the concentration of the CH4 in the atmosphere, defaults to 7.5E-05
    :type ch4conc_atm: float, optional, mol m-3
    :param o2dfc: Diffusion constants (Da) of O2, defaults to 2.02e-05
    :type o2dfc: float, optional, m2 s-1
    :param ch4dfc: Diffusion constants (Da) of CH4, defaults to 1.06e-05
    :type ch4dfc: float, optional, m2 s-1
    :param root_len_den: root length density (Lv, cm root cm-3 soil), defaults to 0.1
    :type root_len_den: float, optional, cm root cm-3 soil
    """
    def __init__(
            self,
            c_aex=24.0,
            eta=400,
            k1=0.33,
            k2=0.44,
            k3=0.22,
            sc_root=3.0e-04,
            o2conc_atm=7.76,
            ch4conc_atm=7.5E-05,
            o2dfc = 2.02e-05,
            ch4dfc = 1.06e-05,
            root_len_den = 0.1, # FIXME: I guess this value.
            ch4sol=1.31, # at 25 °C
            o2sol=1.23,  # at 25 °C
            ke=8.64e+04 # As the time-step in the CERES-Rice model is 1 d, ke takes a value of 8.64 × 104 s.

            ):

        self.c_aex = c_aex 
        self.eta = eta  
        self.k1 = k1 
        self.k2 = k2
        self.k3 = k3 
        self.sc_root = sc_root # the specific conductivity (units: m air (m root)-1) of the root system.
        self.o2conc_atm = o2conc_atm
        self.root_len_den = root_len_den
        self.o2dfc = o2dfc
        self.ch4dfc = ch4dfc
        self.ch4conc_atm = ch4conc_atm
        self.ch4sol = ch4sol
        self.o2sol = o2sol
        self.ke = ke

