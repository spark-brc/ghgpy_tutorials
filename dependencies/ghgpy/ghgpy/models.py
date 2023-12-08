"""
Hartman, M.D., Parton, W.J., Grosso, S.J.D., Easter, M., Hendryx, J., Hilinski, T., Kelly, R., 
Keough, C.A., Killian, K., Lutz, S., Marx, E., McKeown, R., Ogle, S., Ojima, D.S., Paustian, K., 
Swan, A., Williams, S., n.d. 
The Daily Century Ecosystem, Soil Organic Matter, Nutrient Cycling, Nitrogen Trace Gas, 
and Methane Model User Manual, Scientific Basis, and Technical Documentation.

Anaerobic carbohydrate fermentation and methanogenesis occur through the following reactions, 
2(CH2O) → CO2+ CH4 (Conrad 1989, Matthews et al. 2000) or C6H12O6→ 3CO2+ 3CH4 (Huang et al. 1998). 
Both of these equations illustrate that the carbon substrate producing CH4 also produces 
the same carbon equivalent of CO2. 
Hence, a conversion factor C6H12O6_to_CH4 on a mole weight basis of carbohydrate to CH4 is set 
at 0.5 based on the reactions.

"""
import os
import math
import numpy as np
from ghgpy.parms import DCparms, MERESparms, DNDCparms
import pandas as pd


class DCmodel(object):

    def __init__(self, model_dir):
        self.model_dir = model_dir
        # self.beta1, self.Rh, self.fracToExduates, self.sand_cont = self.read_parms()
        self.parms = DCparms()


    def ch4prod(self, sand_cont, eh_t, t_soil, root_c_prod):
        """CH4 production was simulated based on carbon substrate supply and 
        associated influence of Eh and temperature `(Huang et al., 1998)
        <https://doi-org.srv-proxy2.library.tamu.edu/10.1046/j.1365-2486.1998.00129.x>`_.

        .. math::
            \\alpha_1\\times F_{Eh}\\times(C_{soil}+F_T\\times C_{root})

        :param sand_cont: the average sand content fraction (sand, 0.0 - 1.0) in the top 10 cm of soil
        :type sand_cont: float
        :param eh_t: initial redox potential
        :type eh_t: float, mV
        :param t_soil: average soil temperature in the top 10 cm of soil (◦C)
        :type t_soil: float, ◦C
        :param root_c_prod: the previous day's fine root production estimated by
                            the plant production submodel in DayCent
        :type root_c_prod: float, :math:`g\;C\; m^{-2}\cdot d^{-1}`
        :return: CH4Prod is the CH4 production rate.
        :rtype: float, :math:`g\;CH4-C\; m^{-2}\cdot d^{-1}`

        """
        ch4prod_ = (
            self.parms.cvfr_cho_to_ch4 * 
            self.feh(eh_t) * 
            (self.c_soil(sand_cont) + self.f_temp(t_soil) * self.c_root(sand_cont, root_c_prod))
            )
        return ch4prod_

    def c_soil(self, sand_cont):
        """The first step in modeling methanogenesis is to estimate 
        the amount of carbon substrate available for CH4 production. 
        DayCent's methanogenesis submodel includes soil organic matter degradation and 
        rhizodeposition as the sources of carbon.

        :arg sand_cont2: test
        :type sand_cont2: float
        :param sand_cont: _description_
        :type sand_cont: _type_
        :return: _description_
        :rtype: _type_
        """
        c_soil_ = self.parms.cvcf_oc_to_co2 * self.parms.frToCH4 * self.soil_index(sand_cont) * self.parms.hr
        return c_soil_

    def soil_index(self, sand_cont):
        """calculate soil texture index

        Args:
            sand_cont (float): index which is a function of the average sand content fraction 
            (sand, 0.0 - 1.0) in the top 10 cm of soil

        Returns:
            float: soil texture index
        """
        soil_index_ = 0.325+2.25* sand_cont
        return soil_index_
    

    def c_root(self, root_c_prod, sand_cont):
        """The rate of rhizodeposition (Croot, gC m^-2 d^-1) is calculated using 
        the following equation.

        Args:
            frac_to_exudates (float): the fraction of root carbon production contributing to 
            rhizodeposition (range of 0.35-0.6 described in Cao et al. (1995))
            root_c_prod (_type_): the previous day's fine root production estimated by 
            the plant production submodel in DayCent (gC m^-2 d^-1)
        """
        si = self.soil_index(sand_cont)


        c_root_ = self.parms.frac_to_exd * si * root_c_prod
        return c_root_

    def f_temp(self, t_soil):
        """estimate the influence of soil temperature. To simulate CH4 production, 
        we adopted the approach by (Huang et al. 1998) and (Huang et al. 2004).

        :param t_soil: the average soil temperature in the top 10 cm of soil (◦C)
        :type t_soil: float
        :param q10: a temperature coefficient representing the change of a biological or 
                chemical system as a consequence of increasing the temperature by 10°C, 
                and was assumed to be a value of 3.0 (Huang et al. 1998). Defaults to 3.0.
        :type q10: float
        :return: the temperature factor for CH4 production
        :rtype: float
        """
        q10 = self.parms.q10

        if t_soil > 30:
            t_soil = 30
            f_t = q10 **((t_soil-30)/10)
        else:
            f_t = q10 **((t_soil-30)/10)
        return f_t


    def feh(self, eh_t):
        """FEh, a reduction factor for soil redox potential (Eh) (mV), is 
        estimated using the following equations from Huang et al. (1998) and 
        Huang et al. (2004):

        Args:
            eh_t (float, mV): Eht represents the Eh value at time t, and 
            t is the number of days after flooding began or 
            since drainage occurred in the cycle.
        """
        if eh_t >= -150:
            feh_t = math.exp(-1.7*(150+eh_t)/150)
        else:
            eh_t = -150
            feh_t = math.exp(-1.7*(150+eh_t)/150)
        return feh_t


    def get_eh(
            self, sand_cont, waterlevel, eh_init):
        """DEh and AEh (DEH and AEH, fix.100) are differential coefficients that 
        were estimated as 0.16 and 0.23, respectively. 
        The BEhflood (BEHFL, fix.100) is set at a low-limit value of -250 mV, 
        and Behdrain (BEHDR, fix.100) is set to an upper-limit value of 300 mV (Cao et al. 1995, Huang et al. 2004). 
        Soil Eh is a constant value of -20 mV when intermittent irrigation is used in 
        rice paddies as discussed in Huang et al. (2004).

        Args:
            beta1 (float): A fraction (β1) is defined to quantify 
                        the amount of substrate available for methanogens based on 
                        the simulation of heterotrophic respiration by the DayCent model.
                        β1 (CO2CH4, fix.100) is the fraction of Rh converted to CH4 under 
                        anaerobic conditions;
            SI (_type_): _description_
            Rh (float): Rh is heterotrophic respiration from decomposition of organic matter
                        (above- and below-ground structural and metabolic litter and 
                        above- and below-ground SOC pools) (g CO2^-C m^-2 d^-1)
            waterlevel (_type_): _description_
            eh_t (_type_): _description_
            deh (float, optional): _description_. Defaults to 0.16.
            aeh (float, optional): _description_. Defaults to 0.23.
            beh_flood (int, optional): _description_. Defaults to -250.
            beh_drain (int, optional): _description_. Defaults to 300.

        Returns:
            float, mV: current eh
        """


        c_soil = self.c_soil(sand_cont)

        if waterlevel == "flooding":
            eh_now = eh_init - (self.parms.deh * (self.parms.aeh + min(1, c_soil)) * (eh_init - self.parms.beh_flood))
        elif waterlevel == "draining":
            eh_now = eh_init - (self.parms.deh * (self.parms.aeh + 0.7) * (eh_init - self.parms.beh_drain))
        else: # water added via rain and irrigation events
            eh_now = -20
        return eh_now        


    def ch4ep(self, fp, ch4prod):
        """CH4 emission rates through the rice plants (CH4EP) (g CH4-C m^-2d^-1) were simulated

        Args:
            fp (float): fraction of CH4 emitted via rice plants
            ch4prod (float): CH4Prod is total methane production (g CH4-C m^-2d^-1)
        Returns:
            float: transport of CH4 via ebullition to the atmosphere (CH4Ebl)
        """

        ch4ep_ = fp * ch4prod
        return ch4ep_


    def fp(self, aglivc):
        """ get the fraction of CH4 emitted via rice plants

        Args:
            aglivc (float, g C m^-2): the amount of above-ground live C for the crop as simulated by DayCent (g C m^-2)
            tmxbio (float, g biomass m^-2): the maximum aboveground biomass at the end of growing season
            mxch4f (float, optional): MaXimum Fraction of CH4 production emitted by plants. Defaults to 0.55.

        Returns:
            float: the fraction of CH4 emitted via rice plants
        """
        # the multiplier 2.5 (g biomass/g C) converts C to biomass (g biomass m^-2)
        fp_ = self.parms.mxch4f * (1.0 - (aglivc * 2.5/self.parms.tmxbio))**0.25

        return fp_
    
    def ch4ebl(self, tsoil, ch4prod, ch4ep, bglivc):
        """transport of CH4 via ebullition to the atmosphere (CH4Ebl) was also adopted from Huang et al. (2004)

        Args:
            methzr (float): fraction of CH4 emitted via bubbles when there is zero fine root biomass
            tsoil (_type_): 
            ch4prod (_type_): _description_
            ch4ep (_type_): _description_
            mo (float, gC m^-2d^-1): Mo was set to 0.0015 gC m^-2d^-1 (Huang et al. 2004) but is set to 0.002 in DayCent
            mrtblm (float, g biomass m-2): the root biomass that starts to reduce CH4 bubble formation (g biomass m-2)
            bglivc (float, g C m^-2): the amount of fine root C for the crop as simulated by DayCent (g C m^-2)

        Returns:
            float: transport of CH4 via ebullition to the atmosphere (CH4Ebl)
        """
        
        # the multiplier 2.5 (g biomass/g C) converts C to biomass
        # CH4 ebullition is reduced when Tsoil < 2.718282 °C.
        ch4ebl_ = (
            self.parms.frCH4emit_b * 
            (ch4prod - ch4ep - self.parms.mo) * 
            min(np.log(tsoil), 1.0) * 
            (self.parms.mrtblm/(bglivc*2.5))
            )
        return ch4ebl_
    
    def read_inputs(self):
        input_df = pd.read_csv(
            os.path.join(
                self.model_dir, "input_ghg.csv"), na_values=[-9999, ""],
                index_col=0, parse_dates=True)
        return input_df
    

class MERES(object):
    def __init__(self, model_dir):
        self.model_dir =  model_dir
        self.parms = MERESparms()


    def ch4prod(self, pch4prod, o2conc):
        """Actual CH4 production (CH4) in a given soil layer

        :param pch4prod: potential CH4 production
        :type pch4prod: float, :math:`\\frac{mol\;C}{m^3\cdot s}`
        :param o2conc: concentration of O2
        :type o2conc: float, :math:`mol\; m^{-3}`
        :return: Actual CH4 production
        :rtype: float, :math:`mol\; m^{-3}\cdo s^{-1}`
        """
        ch4prod_ = pch4prod / (1 + (self.parms.eta *o2conc))
        return ch4prod_


    def c_root(self, root_wt):
        """the rate of root exudation is calculated as 
        the product of organic compounds per unit of root biomass 
        (depending on the crop growth stage) and 
        the root weight in each soil layer based on results from Lu et al. (1999).

        :param root_wt: root weight
        :type root_wt: float, :math: `kg\;DM\; ha^{-1}`
        :return: rate of root exudation
        :rtype: float, :math:`g\; C\; m^{-2}d^{-1}`
        """

        return 0.4*0.02* root_wt

    def pch4prod(self, aex, subst_c_prod):
        """potential CH4 production (PCH4*)

        :param aex: alternative electron acceptors in oxidized form (:math:`AEX_{ox}`)
        :type aex: float, :math:`mol\; Ceq\; m^{-3}`
        :param subst_c_prod: the rate of substrate-C production
        :type subst_c_prod: float, :math:`mol\; Ceq\; m^{-3}s^{-1}`
        :return: potential CH4 production (PCH4*)
        :rtype: float, :math:`mol\; C\; m^{-3}s^{-1}`
        """

        if aex > self.parms.c_aex:
            pch4prod_ = 0.0
        elif aex > 0.0 and aex < self.parms.c_aex:
            pch4prod_ = min(0.2 * (1-(aex/self.parms.c_aex)), subst_c_prod)
        elif aex == 0:
            pch4prod_ = subst_c_prod
        return pch4prod_

    def ch4oxid(self, pch4prod, ch4conc, o2conc):
        """
        The rate of CH4 consumption (QCH4, mol m-3 s-1) by 
        the methanotrophic bacteria (see equation 2) in 
        a soil layer is given by the Michaelis-Menten equation
        
        .. math:: 
            p=f(x)

        :param pch4prod: potential CH4 production
        :type pch4prod: float, mol Ceq m-3
        :param ch4conc: ch4 concentration
        :type ch4conc: float, mol m-3
        :param o2conc: o2 concentration
        :type o2conc: float, mol m-3
        :return: rate of CH4 consumption
        :rtype: float, mol m-3 s-1

        """
        ch4oxid_ = (
            pch4prod * 
            (ch4conc/(self.parms.k1 + ch4conc)) * 
            (o2conc/(self.parms.k2 + o2conc))
        )
        # self.parms.
        return ch4oxid_
    
    def o2cons(self, o2conc, ch4oxid, pch4prod):
        """Oxygen consumption rate (QO2, mol m-3 s-1)

        :param o2conc: o2 concentration
        :type o2conc: float, mol m-3
        :param ch4oxid: rate of CH4 consumption
        :type ch4oxid: float, mol m-3 s-1
        :param pch4prod: potential CH4 production
        :type pch4prod: float, mol C m-3 s-1
        :return: Oxygen consumption rate
        :rtype: float, mol C m-3 s-1
        """
        o2cons_ = 2*ch4oxid + 2*pch4prod * (o2conc / (self.parms.k3 + o2conc))
        return o2cons_

    # NOTE: Should we distinguish between O2 concentration and its concentration in the gaseous phase?
    def o2flux(self, o2conc):
        """the fluxes of O2

        :param o2conc: o2 concentration
        :type o2conc: float, mol m-3
        :return: the fluxes of O2
        :rtype: float, mol m-3 s-1
        """
        
        o2flux_ = (
            self.parms.sc_root * 
            (self.parms.root_len_den * 10e+04) * 
            self.parms.o2dfc *
            (self.parms.o2conc_atm - o2conc)
        )
        return o2flux_


    def ch4ebl(self, ch4conc_sol):
        """We have modified the algorithm describing ebullition rate from that in the original Arah & Kirk (2000) model by 
        expressing the rate of ebullition (E, mol m-3 s-1) as a function of the concentration of 
        the substance in solution (yw, mol m-3). 
        Currently, there is no temperature dependence of yw* included in the model. 


        :param ch4conc_sol: concentration of the substance in solution
        :type ch4conc_sol: float, mol m-3
        :return: rate of ebullition
        :rtype: float, mol m-3 s-1
        """

        ch4ebl_ = max(
            0, 
            (ch4conc_sol - self.parms.ch4sol)/self.parms.ke
            )
        return ch4ebl_
    

    def o2ebl(self, o2conc_sol):
        """We have modified the algorithm describing ebullition rate from that in the original Arah & Kirk (2000) model by 
        expressing the rate of ebullition (E, mol m-3 s-1) as a function of the concentration of 
        the substance in solution (yw, mol m-3). 
        Currently, there is no temperature dependence of yw* included in the model. 


        :param ch4conc_sol: concentration of the substance in solution
        :type ch4conc_sol: float, mol m-3
        :return: rate of ebullition
        :rtype: float, mol m-3 s-1
        """

        o2ebl_ = max(
            0, 
            (o2conc_sol - self.parms.o2sol)/self.parms.ke
            )
        return o2ebl_
    


class DNDC(object):
    def __init__(self, model_dir):
        self.model_dir =  model_dir
        self.parms = DNDCparms()

    def ch4prod(self, ava_c, ft_temp):
        """
        .. math::
            CH_{4}p = a * AC * Ft

        :param ava_c: the available carbon concentration
        :type ava_c: float, :math:`kg C/ha`
        :param ft_temp: the temperature factor,
        :type ft_temp: float
        :return: CH4 production rate
        :rtype: float, :math:`\\frac {kg C}{ha*d}`
        """

        ch4prod_ = self.parms.a * ava_c * ft_temp
        return ch4prod_
    
    def ft_temp(self, t_soil):
        """the temperature factor
        .. math::
            Ft = b*e^{0.2424*T}

        :param t_soil: the temperature (°C)
        :type t_soil: float, °C
        :const b: constant parameter
        :type float
        :return: temperature factor
        :rtype: float
        """
        ft_temp_ = self.parms.b * math.exp(0.2424*t_soil)
        return ft_temp_

    def ch4oxid(self, ch4conc, eh):
        ch4oxid_ = ch4conc * math.exp(8.6711*eh/1000)
        return ch4oxid_
    
    def ch4ep(self, ch4prod, aere):
        ch4ep_ = 0.5 * ch4prod * aere
        return ch4ep_
    
    def aere(self, pgi):
        """ Calculate AERE (Aerobic Respiration Enhancement) based on the Polynomial Equation.

        :param pgi: Plant Growth Index.
        :type pgi: float
        :return: The calculated AERE value.
        :rtype: float
        """
        # Define coefficients for the polynomial equation
        a = -0.0009
        b = 0.0047
        c = -0.883
        d = 1.9863
        e = -0.3795
        f = 0.0251

        # Calculate AERE using the polynomial equation
        aere_ = a * pgi**5 + b * pgi**4 + c * pgi**3 + d * pgi**2 + e * pgi + f
        return aere_
    
    def pgi(self, stdate, eddate, curdate):
        """plant growth index

        :param stdate: start date
        :type stdate: str, dateformat
        :param dsp: days since planting
        :type dsp: int
        :param sds: season days
        :type sds: int
        :return: plant growth index
        :rtype: float
        """
        stdate =  pd.to_datetime(stdate)
        eddate = pd.to_datetime(eddate)
        dsp = abs((curdate - stdate).days)
        sds = abs((eddate - stdate).days) 

        pgi_ = dsp / sds
        return pgi_
    
    def ch4ebl(self, ch4prod, poro, ft_temp2, aere):
        """Calculate CH4 emission via ebullition based on various factors.

        :param ch4prod: Total methane production.
        :type ch4prod: float
        :param poro: Soil porosity.
        :type poro: float
        :param ft_temp2: Temperature-based factor.
        :type ft_temp2: float
        :param aere: aerenchyma factor
        :type aere: float
        :return: CH4 emission via ebullition.
        :rtype: float
        """

        ch4ebl_ = 0.25 * ch4prod* poro * ft_temp2 * (1 - aere)
        return ch4ebl_

    def ft_temp2(self, t_soil):
        """Calculate a temperature-based factor (ft_temp2) based on a polynomial equation.

        :param t_soil: The average soil temperature in °C.
        :type t_soil: float
        :return: The calculated temperature-based factor (ft_temp2).
        :rtype: float
        """

        # Coefficients for the polynomial equation
        a = -0.1687
        b = 1.167
        c = -2.0303
        d = 1.042

        # Calculate ft_temp2 using the polynomial equation
        t_soil_normalized = 0.1 * t_soil
        ft_temp2_ = a * t_soil_normalized**3 + b * t_soil_normalized**2 + c * t_soil_normalized + d

        return ft_temp2_
    
    def dffs_rate(self, ch4prod, poro, t_soil):
        dffs_rate_ = 0.01 * ch4prod * t_soil * poro
        return dffs_rate_
    

    def eh(self, e0, t_soil, num_elect):
        eh_ = e0 + (self.parms.rconst * t_soil / (num_elect * self.parms.fconst))
        return eh_

    def read_inputs(self):
        input_df = pd.read_csv(
            os.path.join(
                self.model_dir, "input_ghg.csv"), na_values=[-9999, ""],
                index_col=0, parse_dates=True)
        return input_df
    




class DNDCdaily(object):
    def __init__(self, model_dir):
        self.model_dir =  model_dir
        self.parms = DNDCparms()

    def ftm(self, t_soil):
        ftm_ = math.exp(
            0.33* (t_soil - 23) / 
            (1 + math.exp(0.33*(t_soil-23)))
            )
        return ftm_
    
    def fphm(self, ph):
        """Calculate a pH-based factor for methane production. 
        This function computes a factor that is based on pH and 
        can be used in methane production models.

        :param ph: The pH value
        :type ph: float
        :return: The calculated pH-based factor.
        :rtype: float
        """
        numerator = (ph - 5.5) * (ph - 9.0)
        denominator = (ph - 5.5) * (ph - 9.0) - (ph - 7.5) ** 2

        # Avoid division by zero and return 0 if the denominator is zero
        fphm_ = numerator / denominator if denominator != 0 else 0.0

        return fphm_
    
    def feh(self, eh):
        if eh <= -200:
            feh_ = 1
        else:
            feh_ = 0
        return feh_
    
    def c_pool(self, c_sol, bm_root):
        """Calculate the carbon pool for methane production (kg C ha-1).
        This function computes the carbon pool available for methane production based on 
        soluble carbon and root biomass.

        :param c_sol: The soluble carbon content in the layer (kg C ha-1).
        :type c_sol: float, kg C ha-1
        :param bm_root: the root biomass (g m-2)
        :type bm_root: float, g m-2)
        :return: The calculated carbon pool for methane production (kg C ha-1).
        :rtype: kg C ha-1
        """
        # Constants from parameters (consider defining them as class-level constants)
        fr_rsi = self.parms.fr_rsi
        fr_root = self.parms.fr_root

        # Calculate the carbon pool
        c_pool_ = c_sol + fr_rsi * fr_root * bm_root * 4

        return c_pool_
    
    def ch4prod(self, c_pool, ftm, feh, fphm):
        ch4prod_ = 0.47 * c_pool * ftm * feh * fphm
        return ch4prod_

    def ch4oxid(self, ch4prod, aere):
        ch4oxid_ = ch4prod * (0.5 + 0.5 * aere)
        return ch4oxid_

    def aere(self, bm_root):
        aere_ =  bm_root / 1000
        return aere_

    def ch4e(self, ch4prod, ch4oxid):
        ch4e_ = ch4prod - ch4oxid
        return ch4e_
    
    def read_inputs(self):
        input_df = pd.read_csv(
            os.path.join(
                self.model_dir, "input_ghg.csv"), na_values=[-9999, ""],
                index_col=0, parse_dates=True)
        return input_df