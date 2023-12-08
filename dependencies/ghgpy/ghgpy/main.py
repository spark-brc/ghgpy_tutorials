# -*- coding: utf-8 -*-
"""
.. math::
    \\frac{ \sum_{t=0}^{N}f(t,k) }{N}

"""

import os
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(".."))

from ghgpy.models import DCmodel
from ghgpy.utils import ObjFns
from ghgpy.utils import PostPr
import pandas as pd
import os
from ghgpy.analyzer import (
    plot_oom, plot_tseries_ch4, plot_tseries_ch4_tot, plot_violin, plot_oot)
from ghgpy.runs import run_dc_multi, run_dndc_daily, run_dndc



def run_dndc_daily_model(wd):
    run_dndc_daily(wd)
    so_df = PostPr(wd).get_ch4_so_df(outfnam="ch4_multi_dndc_daily.out")
    # viz
    plot_tseries_ch4(so_df, simnam="ch4es", height=3, dot=False)

def run_dndc_model(wd):
    run_dndc(wd)
    so_df = PostPr(wd).get_ch4_so_df(outfnam="ch4_multi_dndc.out")
    # viz
    plot_tseries_ch4(so_df, simnam="ch4e_tot", height=3, dot=False)


def run_dc_model(wd):
    run_dc_multi(wd)
    so_df = PostPr(wd).get_ch4_so_df(outfnam="ch4_multi_dc.out")
    plot_tseries_ch4(so_df, height=3, dot=False, fignam="dc_out.png")

def run_md_viz(wd):
    # run model
    run_dc_multi(wd)
    # post-processing
    m1 = PostPr(wd)


    so_df = m1.get_ch4_so_df(outfnam="ch4_multi_dndc.out")
    # so_df = m1.get_ch4_so_df(outfnam="ch4_multi.out")
    print(so_df)

    # viz
    plot_tseries_ch4(so_df, dot=False)
    


def viz_tot(wd):
    df = pd.read_excel(
        os.path.join(wd, "tot_3m.xlsx"),
        # index_col=0,
        # parse_dates=True,
        comment="#")
    df = df.loc[df["cont"]=="flooding"]
    dff = df.loc[:, ["ch4_obd", "dc", "dndc", "dndc_d"]]
    # plot_violin(wd, dff)
    plot_oot(wd, dff)



if __name__ == "__main__":
    wd = "D:\\Projects\\Tools\\ghgpy\\models"
    os.chdir(wd)
    run_dc_model(wd)

