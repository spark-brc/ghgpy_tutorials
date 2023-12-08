from ghgpy.models import DCmodel, DNDCdaily, DNDC
import os
import pandas as pd


SFMT = lambda x: "{0:<10s} ".format(str(x))
IFMT = lambda x: "{0:<10d} ".format(int(x))
FFMT = lambda x: "{0:<15.10E} ".format(float(x))


def run_dc_single(wd):
    m1 = DCmodel(wd)
    sand_cont = 0.038
    root_c_prod = 6

    indf = m1.read_inputs()
    indf.dropna(axis=0, inplace=True)

    outdf = pd.DataFrame(
        columns=['ch4prod','ch4ep','ch4ebl'],
        index=indf.index)
    outdf.index.name = "date"
    # outdf.index = outdf.index.strftime("%Y-%m-%d")

    for i in indf.index:
        print(indf.loc[i])
        ch4prod = m1.ch4prod(
            sand_cont, float(indf.loc[i, "eh"]), 
            float(indf.loc[i, "tsoil"]), root_c_prod)
        outdf.loc[i, "ch4prod"] = ch4prod
    outdf.to_csv(os.path.join(wd,"ch4_output.csv"))


def run_dc_multi(wd):
    m1 = DCmodel(wd)
    # read inputs
    sand_cont = 0.038
    root_c_prod = 80
    aglivc = 100 # the amount of above-ground live C for the crop as simulated by DayCent (g C m−2)
    bglivc = 100 # the amount of fine root C for the crop as simulated by DayCent (g C m−2)
    #---------------
    indf = m1.read_inputs()
    indf.dropna(axis=0, inplace=True)
    conts = indf["condition"].unique()
    dff = pd.DataFrame()

    for cont in conts:
        contdf = indf.loc[indf["condition"]==cont]
        ch4prods = []
        ch4eps = []
        ch4ebls = []
        for i in contdf.index:
            # print(contdf.loc[i])
            # calculate CH4 production
            ch4prod = m1.ch4prod(
                sand_cont, float(contdf.loc[i, "eh"]), 
                float(contdf.loc[i, "tsoil"]), root_c_prod)
            ch4prods.append(ch4prod)
            # calculate CH4 emission by plant
            fp = m1.fp(aglivc)
            ch4ep = m1.ch4ep(fp, ch4prod)
            ch4eps.append(ch4ep)
            # calculate CH4 emission by ebullision
            ch4ebl = m1.ch4ebl(
                float(contdf.loc[i, "tsoil"]), ch4prod, ch4ep, bglivc)
            ch4ebls.append(ch4ebl)

        getdf = pd.DataFrame(
            {"ch4prods":ch4prods, "ch4eps":ch4eps, "ch4ebls":ch4ebls}, 
            index=contdf.index)
        getdf["ch4e_tot"] = getdf["ch4eps"] + getdf["ch4ebls"]
        getdf['cont'] = cont
        getdf["ch4_obd"] = contdf.loc[:, "ch4_obd"]
        dff = pd.concat([dff, getdf], axis=0)
    dff.insert(0, "date", dff.index.date)         
    with open(os.path.join(wd, "ch4_multi_dc.out"), "w") as f:
        fmts = [SFMT, FFMT, FFMT, FFMT, FFMT, FFMT, SFMT] 
        f.write("# created by ghgpy\n")
        f.write(dff.loc[:, [
            "date", "ch4_obd", "ch4prods", "ch4eps", "ch4ebls", "ch4e_tot", "cont"]].to_string(
                                                        col_space=0,
                                                        formatters=fmts,
                                                        index=False,
                                                        header=True,
                                                        justify="left"))

    # outdf.index = outdf.index.strftime("%Y-%m-%d")
    print("run complete ...")


def run_dndc_daily(wd):
    m1 = DNDCdaily(wd)
    ph = 6.72
    root_biomass =2
    c_sol = 2
    # read inputs
    indf = m1.read_inputs()
    indf.dropna(axis=0, inplace=True)
    conts = indf["condition"].unique()
    dff = pd.DataFrame()

    for cont in conts:
        contdf = indf.loc[indf["condition"]==cont]
        c_pools = []
        fehs = []
        ch4prods = []
        ch4oxids = []
        ch4es = []
        for i in contdf.index:
            # print(contdf.loc[i])
            # calculate c_pool
            c_pool = m1.c_pool(c_sol, root_biomass)
            c_pools.append(c_pool)
            # calculate feh
            feh = m1.feh(float(contdf.loc[i, "eh"]))
            fehs.append(feh)
            # calculate ftm
            ftm = m1.ftm(float(contdf.loc[i, "tsoil"]))
            # calculate fphm
            fphm = m1.fphm(ph)
            # calculate ch4prod
            ch4prod = m1.ch4prod(c_pool, ftm, feh, fphm)
            ch4prods.append(ch4prod)
            # calculate aere
            aere = m1.aere(root_biomass)
            # calcualte ch4oxid
            ch4oxid = m1.ch4oxid(ch4prod, aere)
            ch4oxids.append(ch4oxid)
            # calcualte emission
            ch4e = m1.ch4e(ch4prod, ch4oxid)
            ch4es.append(ch4e)
        getdf = pd.DataFrame(
            {"ch4prods":ch4prods, "ch4oxids":ch4oxids, "ch4es":ch4es}, 
            index=contdf.index)
        # getdf["ch4e_tot"] = getdf["ch4eps"] + getdf["ch4ebls"]
        getdf['cont'] = cont
        getdf["ch4_obd"] = contdf.loc[:, "ch4_obd"]
        dff = pd.concat([dff, getdf], axis=0)
    dff.insert(0, "date", dff.index.date)
    with open(os.path.join(wd, "ch4_multi_dndc_daily.out"), "w") as f:
        fmts = [SFMT, FFMT, FFMT, FFMT, FFMT, SFMT] 
        f.write("# created by ghgpy\n")
        f.write(dff.loc[:, [
            "date", "ch4_obd", "ch4prods", "ch4oxids", "ch4es", "cont"]].to_string(
                                                        col_space=0,
                                                        formatters=fmts,
                                                        index=False,
                                                        header=True,
                                                        justify="left"))
    print("run complete ...")

def run_dndc(wd):
    m1 = DNDC(wd)
    ava_c = 0.07
    # pg_idx = 0.5
    poro = 0.01
    stdate = "5/19/2022"
    eddate = "10/5/2022"
    # root_biomass =1
    # c_sol = 1
    # read inputs
    indf = m1.read_inputs()
    indf.dropna(axis=0, inplace=True)
    conts = indf["condition"].unique()
    dff = pd.DataFrame()

    for cont in conts:
        contdf = indf.loc[indf["condition"]==cont]
        # c_pools = []
        # fehs = []
        ch4prods = []
        ch4oxids = []
        ch4eps = []
        ch4ebls = []
        for i in contdf.index:
            tempf01 = m1.ft_temp(float(contdf.loc[i, "tsoil"]))
            ch4prod = m1.ch4prod(ava_c, tempf01) * 0.1
            ch4prods.append(ch4prod)
            # NOTE: use ch4p as ch4 concentration?
            ch4oxid = m1.ch4oxid(ch4prod, float(contdf.loc[i, "eh"]))
            ch4oxids.append(ch4oxid)
            pg_idx = m1.pgi(stdate, eddate, i)
            aere = m1.aere(pg_idx)
            ch4ep = m1.ch4ep(ch4prod, aere)
            ch4eps.append(ch4ep)
            tempf02 = m1.ft_temp2(float(contdf.loc[i, "tsoil"]))
            ch4ebl = m1.ch4ebl(ch4prod, poro, tempf02, aere)
            ch4ebls.append(ch4ebl)

        getdf = pd.DataFrame(
            {"ch4prods":ch4prods, "ch4eps":ch4eps, "ch4ebls":ch4ebls}, 
            index=contdf.index)
        getdf["ch4e_tot"] = getdf["ch4eps"] + getdf["ch4ebls"]
        getdf['cont'] = cont
        getdf["ch4_obd"] = contdf.loc[:, "ch4_obd"]
        dff = pd.concat([dff, getdf], axis=0)
    dff.insert(0, "date", dff.index.date)
    with open(os.path.join(wd, "ch4_multi_dndc.out"), "w") as f:
        fmts = [SFMT, FFMT, FFMT, FFMT, FFMT, FFMT, SFMT]
        f.write("# created by ghgpy\n")
        f.write(dff.loc[:, [
            "date", "ch4_obd", "ch4prods", "ch4eps", "ch4ebls", "ch4e_tot", "cont"]].to_string(
                                                        col_space=0,
                                                        formatters=fmts,
                                                        index=False,
                                                        header=True,
                                                        justify="left"))
    print("run complete ...")
