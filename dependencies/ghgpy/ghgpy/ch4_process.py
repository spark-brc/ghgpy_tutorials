# -*- coding: utf-8 -*-
class Decomposition(object):
    """
    FOM (crop residues, green manure, and more stable organic or humic pool (HUM))
    The CERES-Rice model requires input data (FOM): 
    .. Amount of straw added,
    .. C/N
    .. Depth of incorporation
    .. root residue from previous crop

    Initialization of the HUM pool is calculated from the soil organic carbon as specified in the soil data file.
    """

    def __init__(self):
        """
        MERES Decay rate (Seligman and Keulen, 1981)
        daily decay rate
        """
        self.dr_carb = 0.2 # carbohydrate
        self.dr_cel = 0.05 # cellulose
        self.dr_lig = 0.0095 # lignin

        """
        DayCent Decomposition rate
        per week
        Ki: state variables
        """
        self.k01 = 0.076 # structural soil surface litter
        self.k02 = 0.28 # metaboilc soil surface litter
        self.k03 = 0.094 # structural soil litter
        self.k04 = 0.35 # metaboilc soil litter
        self.k05 = 0.14 # active soil fraction
        self.k06 = 0.0038 # slow soil fraction
        self.k07 = 0.00013 # passive soil fraction

    def actual_dr(self, Qp): #dQp/dt, kg C ha^-1 d^-1
        """
        Eq:
            dQp/dt = QpRp(max) * f(Ts) * g(theta_s) * h(k_s)
        .. Qp (kg C ha-1): amount of organic matter in pool p on the day
        .. Rp(max): decay constants
        .. f(Ts): dimesionless multiplier for soil temperature (C)
        .. g(theta_s): dimesionless multiplier for soil moisture (m3 water m^-3 soil)
        .. h(k_s): pool C/N (kg C kg^-1 N)
        Returns:
            _type_: float
            The total amount of carbon released by decay of organic matter on a given day
        """
        Rps = [self.dr_carb, self.dr_cel, self.dr_lig]
        tot_Op = 0
        for Rp in Rps:  
            dOp = Qp * Rp * self.dr_temp() * self.dr_moisture() * self.dr_cnpool()
            tot_Op += dOp
        return tot_Op
    
    def decomp_rate(self):
        """
        DayCent model
        eq:
            \sfrac{{dC}_i}{d_t}=K_i * M_d * T_d * C_i
        """

    def dr_temp(self):
        return 
    
    def dr_moisture(self):
        return 
    
    def dr_cnpool(self):
        return

    def rel_Cdecay(self):
        return


class rhizoDeposit(object):
    """
    The contribution to CH4 production of organic matter originating from 
    living rice plants through root exudates and root death, 
    collectively referred to as rhizodeposition.

    - The peak in emission rates commonly observed toward the end of the growing season was ascribed by 
    these authors to be due to the increase in decaying root tissue or root exudates after flowering.


    Args:
        object (_type_): _description_
    """
    def __init__(self):
        """

        beta_2 is the fraction of root carbon production contributing to rhizodeposition, 
        which was set to an intermediate value of 0.45 from a range of 0.35–0.6 described in Cao et al. (1995).
        Args:
            RC (_type_): _description_
        """
        self.beta2 = 0.45


    def root_exudation(self, w_root):
        """
        In MERES, the rate of root exudation (g\ C\ m^{-2}d^{-1}) is calculated as 
        the product of organic compounds per unit of root biomass (depending on the crop growth stage) and 
        the root weight in each soil layer based on results from Lu et al. (1999).        

        Args:
            w_root (float, Kg DM ha^-1): existing root dry weight
            Root death is assumed to be a constant 2% of existing root dry weight 

        Returns:
            float: the rate of root exudation (g C m^{-2}d^{-1})
        """
        
        rCroots = 0.4 * 0.02 * w_root
        return rCroots
    
    def root_exud_dc(self, RC):
        """_summary_

        Args:
            RC (_type_): _description_
        """
        Croot = self.beta2 * RC


class electAcceptor(object):
    """
    Unless Oxygen (O2) is depleted in the soil, O2 acts as the sole electron acceptor for microbial respiration. 
    After rice soil is flooded, dissolved O2 in the floodwater and soil is consumed rapidly. 
    The need for electron acceptors by anaerobic organisms results in 
    the reduction of a number of other oxidized species of ions in the soil. 
    Reductions of NO3- to NO2-, N2O to N2, Mn4+ to Mn2+, Fe3+ to Fe2+, SO42- to S2- all resulting in 
    the production of CO2, and finally CO2 to CH4, occur sequentially, 
    provided available carbon sources exist (Patrick and Delaune, 1977). 
    CH4 production will occur after most of the soil's NO3-, N2O, Mn4+, Fe3+, and 
    SO42- ions have been reduced.

    Args:
        object (_type_): _description_
    """

    def __init__(self):
        pass

    def aex_ox_meres(self):
        """
        The MERES model assumed the presence in the soil of a pool of alternative electron acceptors in 
        oxidized form ({AEX}_{ox}), which reacts with the substrate C from decomposition to form CO2, 
        becoming reduced in the process ({AEX}_{red}). 
        It specifies the quantity of AEA present in moles of C equivalents m-3, 
        assuming a 1:1 stoichiometric relationship between substrate C and the AEA, i.e.,
        \left({CH}_2O\right)+{AEX}_{ox}\rightarrow{CO}_2+{AEX}_{red}+{2H}^+
        - 
        """
    
    def pCH4_aex(self, aex_ox, rSC, aea_ox_max=24.0):
        if aex_ox > aea_ox_max:
            pCH4 = 0 
        elif aea_ox_max > aex_ox:
            pCH4 = min([0.2*(1-aex_ox/aea_ox_max), rSC])


class Oxidation(object):
    """
    MERES
    - uses the Michaelis-Menten equation to simulate the rate of CH4 oxidation (QCH4, mol m-3 s-1) by the methanotrophic bacteria.
    - Eq: Q_{CH4}=P_{CH4} * \ast\bullet[CH4](k1+[CH4])∙[O2](k1+[O2])
    .. P_ch4a: potential rate of methanogenesis defined previously, mol/(m^3*s)
    .. bracket_ch4: concentration of CH4, mol/m^3
    .. bracket_o2: concentration of O2, mol/m3
    .. k1: 0.33, Michaelis-Menten constants (units: mol m-3) for a dual-substrate reaction.
    .. k2: 0.44, Michaelis-Menten constants (units: mol m-3) for a dual-substrate reaction.
    .. k3: 0.22, Michaelis-Menten constants (units: mol m-3) for a dual-substrate reaction.
    """

    def __init__(self):
        self.k1 = 0.33
        self.k2 = 0.44
        self.k3 = 0.22

    def ox_ch4(self, p_ch4, conc_ch4, conc_o2):
        q_ch4 = p_ch4 * (conc_ch4/(self.k1+conc_ch4) * conc_o2/(self.k2 + conc_o2))
        return q_ch4 


class Ebullition(object):
    """
    MERES
    .. E: rate of ebullition, mol/(m^3*s)
    .. y_w: concentration of the substance in solution, mol/m^3
    .. y_wa: solubility of the substance in water, mol/m^3
            Values of y_wa (at 25 °C) used were 1.23 mol m^-3 and 
            1.31 mol m^-3 for O2 and CH4, respectively.      
    .. k_e: a constant (units: s) equal to the timestep of the simulation.
            As the time-step in the SWAT+ model is 1 d, 
            k_e takes a value of 86400s.
    - Currently, there is no temperature dependence of y_wa included in the model, 
        but this could be incorporated in future versions.
    
    """

    def __init__(self, y_wa=1.23, k_e=86400):
        self.y_wa = y_wa
        self.k_e = k_e

    def e_ch4(self, y_w):
        e_ch4 = max([0, (y_w - self.y_wa)/self.k_e])
        return e_ch4


class Aerenchyma(object):
    """
    MERES
    - assumes that the conductance of the plant pathway is proportional to the root length density 
    (Lv, cm root cm-3 soil) present in each soil layer, 
    such that λ= λr Lv, where λr represents the specific conductivity (units: m air (m root)-1) of 
    the root system. 
    Thus, the flux (F, mol m-3 s-1) for each substance (O2 or CH4) is given by    
    - Eq: F=\lambda_r(L_v * {10}^4)D_a(y_{a0}-y_a)
    .. λ: the conductance of the pathway through the plant, m air m-3 soil
    .. λr: specific conductivity (units: m air (m root)-1) of the root system.
    .. L_v: the root length density (Lv, cm root cm-3 soil) present in each soil layer
    .. D_a: the diffusivity of the respective substance through air, m2 s^-1
    .. y_a0: concentration of the respective substance in 
            the atmosphere (O2: 7.76 mol m^-3; CH4: 7.5 * 10-5 mol m^-3),
    .. y_a: concentration in the gaseous phase in each soil layer.
    
    """

    def __init__(self):
        self.y_a0_o2 = 7.76
        self.y_a0_ch4 = 7.5 * 0.00001
        self.D_a_o2 = 2.02 * 0.00001
        self.D_a_ch4 = 1.06 * 0.00001

    def a_ch4(self, y_a, l_v):
        a_ch4 = y_a * (l_v * 10000) * self.D_a_ch4 * (self.y_a0_ch4 - y_a)
        return a_ch4


def c_soil(conv_fac, SI, beta_1, Rh):
    """The first step in modeling methanogenesis is to estimate 
    the amount of carbon substrate available for CH4 production. 
    DayCent's methanogenesis submodel includes soil organic matter degradation and 
    rhizodeposition as the sources of carbon.

    Args:
        conv_fac (float): _description_
        SI (_type_): _description_
        beta_1 (float): _description_
        Rh (_type_): _description_

    Returns:
        _type_: _description_
    """
    c_soil = conv_fac * beta_1 * SI * Rh
    return c_soil

def SI(send_cont_frac):
    """calculate soil texture index

    Args:
        send_cont_frac (float): index which is a function of the average sand content fraction 
        (sand, 0.0 - 1.0) in the top 10 cm of soil

    Returns:
        float: soil texture index
    """
    SI = 0.325+2.25*send_cont_frac
    return SI





def c_root(conv_fac_exudate, rootCprod, send_cont_frac):
    conv_fac_exudate * SI(send_cont_frac) * rootCprod