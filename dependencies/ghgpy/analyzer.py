import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from ghgpy.utils import ObjFns
import numpy as np
import pandas as pd
import os
from ghgpy.models import DCmodel


def plot_oo(df, simnam=None,numcols=1, fsize=8):
    m1 = ObjFns()
    if simnam is None:
        simnam = "ch4e_tot"
    fig, ax = plt.subplots(figsize=(6,5))
    # colors = cm.tab20(np.linspace(0, 1, len(df.site_name.unique())))
    fmax = df.loc[:, ["ch4prods", "ch4_obd"]].max().max()
    fmin = df.loc[:, ["ch4prods", "ch4_obd"]].min().min()
    x_val = df.loc[:, "ch4prods"].tolist()
    y_val = df.loc[:, "ch4_obd"].tolist()
    correlation_matrix = np.corrcoef(x_val, y_val)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    m, b = np.polyfit(x_val, y_val, 1)
    ax.plot(np.array(x_val), (m*np.array(x_val)) + b, 'k', label='_nolegend_')
    ax.text(
            0.05, 0.9,
            f'$R^2:$ {r_squared:.3f}',
            horizontalalignment='left',
            bbox=dict(facecolor='gray', alpha=0.2),
            transform=ax.transAxes
            )
    ax.text(
            0.95, 0.05,
            f'$y={m:.2f}x{b:.2f}$',
            horizontalalignment='right',
            # bbox=dict(facecolor='gray', alpha=0.2),
            transform=ax.transAxes
            )
    ax.scatter(df[simnam], df.ch4_obd,  alpha=0.7)
    rsq_val = round(m1.rsq(df[simnam], df.ch4_obd), 3)
    rmse_val = round(m1.rmse(df[simnam].values, df.ch4_obd.values), 3)
    pbias_val = round(m1.pbias(df[simnam].values, df.ch4_obd.values), 3)

    # lgds = []
    # for tn, c in zip(df.site_name.unique(), colors):
    #     sdf = df.loc[df['site_name'] == tn]
    #     ax.scatter(
    #         sdf.somsc_sim, sdf.obd, 
    #         color = c, 
    #         alpha=0.7)
    #     rsq_val = round(rsquared(sdf.obd, sdf.somsc_sim), 3)
    #     rmse_val = round(rmse(sdf.obd, sdf.somsc_sim), 3)
    #     lgds.append(f"{tn} (rsq:{rsq_val}, rmse:{rmse_val})")
    ax.plot([fmin, fmax], [fmin, fmax], 'k--', alpha=0.2)
    ax.set_xlabel("Simulated CH4 ($CH4-C$ $m^{-2} d^{-1}$)")
    ax.set_ylabel("Observed CH4 ($CH4-C$ $m^{-2} d^{-1}$)")
    # plt.legend(
    #     lgds, 
    #     bbox_to_anchor=(1.05, 1.05), ncols=numcols, fontsize=fsize)
    # fig.tight_layout()
    plt.savefig("plot_oneToOne.jpg", dpi=300, bbox_inches="tight")
    plt.show()


def plot_oom(df, simnam=None, obdnam=None, numcols=1, fsize=8):
    # plot one to one multi sites or scenarios
    m1 = ObjFns()
    if simnam is None:
        simnam = "ch4e_tot"
    if obdnam is None:
        obdnam = "ch4_obd"
    fig, ax = plt.subplots(figsize=(6,5))
    colors = cm.rainbow(np.linspace(0, 1, len(df.cont.unique())))
    fmax = df.loc[:, [simnam, obdnam]].max().max()
    fmin = df.loc[:, [simnam, obdnam]].min().min()
    x_val = df.loc[:, simnam].tolist()
    y_val = df.loc[:, obdnam].tolist()
    correlation_matrix = np.corrcoef(x_val, y_val)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2

    m, b = np.polyfit(x_val, y_val, 1)
    ax.plot(np.array(x_val), (m*np.array(x_val)) + b, 'k', label='_nolegend_')
    
    rsq_val = round(m1.rsq(df[simnam], df.ch4_obd), 3)
    rmse_val = round(m1.rmse(df[simnam].values, df.ch4_obd.values), 3)
    pbias_val = round(m1.pbias(df[simnam].values, df.ch4_obd.values), 3)
    ax.text(
            0.05, 0.9,
            f'$R^2:$ {r_squared:.3f}',
            horizontalalignment='left',
            bbox=dict(facecolor='gray', alpha=0.2),
            transform=ax.transAxes
            )
    ax.text(
            0.3, 0.9,
            f'$R^2:$2 {rsq_val:.3f}',
            horizontalalignment='left',
            bbox=dict(facecolor='gray', alpha=0.2),
            transform=ax.transAxes
            )
    ax.text(
            0.7, 0.9,
            f'$rmse:$ {rmse_val:.3f}',
            horizontalalignment='left',
            bbox=dict(facecolor='gray', alpha=0.2),
            transform=ax.transAxes
            )
    ax.text(
            0.95, 0.05,
            f'$y={m:.2f}x{b:.2f}$',
            horizontalalignment='right',
            # bbox=dict(facecolor='gray', alpha=0.2),
            transform=ax.transAxes
            )
    # ax.scatter(df[simnam], df.ch4_obd,  alpha=0.7)
    lgds = []
    for tn, c in zip(df.cont.unique(), colors):
        sdf = df.loc[df['cont'] == tn]
        ax.scatter(
            sdf[simnam], sdf.ch4_obd, 
            color = c, 
            alpha=0.7)
        rsq_val = round(m1.rsq(sdf[simnam], sdf.ch4_obd), 3)
        rmse_val = round(m1.rmse(sdf[simnam].values, sdf.ch4_obd.values), 3)
        lgds.append(f"{tn} (rsq:{rsq_val}, rmse:{rmse_val})")
    ax.plot([fmin, fmax], [fmin, fmax], 'k--', alpha=0.2)
    ax.set_xlabel("Simulated CH4 ($CH4-C$ $m^{-2} d^{-1}$)")
    ax.set_ylabel("Observed CH4 ($CH4-C$ $m^{-2} d^{-1}$)")
    plt.legend(
        lgds, 
        # bbox_to_anchor=(1.05, 1.1),
        ncols=numcols, fontsize=fsize
        )
    fig.tight_layout()
    plt.savefig("plot_oom.jpg", dpi=300, bbox_inches="tight")
    plt.show()


def plot_tseries_ch4(
                    df, simnam=None, obdnam=None, width=10, height=4, dot=True,
                    fignam=None
                    ):
    if simnam is None:
        simnam = "ch4e_tot"
    if obdnam is None:
        obdnam = "ch4_obd"
    if fignam is None:
        fignam = "sim_obd_ch4.png"


    ogs = df.cont.unique()
    fig,axes = plt.subplots(len(ogs),1,figsize=(width,height*len(ogs)))
    ogs.sort()
    # for each observation group (i.e. timeseries)
    for ax,og in zip(axes,ogs):
        # get values for x axis
        oobs = df.loc[df.cont==og,:].copy()
        tvals = oobs.date.values
        # onames = oobs.obsnme.values
        if dot is True:
            # plot prior
            ax.scatter(tvals, oobs.loc[:, simnam],color="gray",s=30, alpha=0.5)
            # plot posterior
            ax.scatter(tvals, oobs.loc[:, obdnam],color='red',s=30).set_facecolor("none")
        if dot is False:
            # plot prior
            ax.plot(tvals, oobs.loc[:, simnam],color="gray", alpha=0.5)
            # plot posterior
            ax.scatter(tvals, oobs.loc[:, obdnam],color='red',s=30).set_facecolor("none")
            # plot measured+noise 
        ax.tick_params(axis='x', labelrotation=45)
        ax.margins(x=0.01)
        ax.set_title(og,loc="left")
    fig.tight_layout()
    plt.savefig(fignam, dpi=300, bbox_inches="tight")
    plt.show()


def plot_tseries_ch4_tot(
                    df, obdnam=None, width=10, height=4, dot=False,
                    fignam=None
                    ):
    if obdnam is None:
        obdnam = "ch4_obd"
    if fignam is None:
        fignam = "sim_obd_ch4_tot.png"
    ogs = df.cont.unique()
    fig,axes = plt.subplots(len(ogs),1,figsize=(width,height*len(ogs)), sharex=True)
    ogs.sort()
    # for each observation group (i.e. timeseries)
    for ax,og in zip(axes,ogs):
        # get values for x axis
        oobs = df.loc[df.cont==og,:].copy()
        tvals = oobs.date.values
        # tvals = pd.to_datetime(oobs['date']).dt.strftime("%Y-%m-%d")

        # onames = oobs.obsnme.values
        if dot is False:
            # plot prior
            ax.plot(tvals, oobs.loc[:, ["dc", "dndc", "dndc_d"]], alpha=0.5)
            # plot posterior
            ax.scatter(tvals, oobs.loc[:, obdnam],color='red',s=30).set_facecolor("none")
            # plot measured+noise 
        ax.tick_params(axis='x', labelrotation=45)
        ax.margins(x=0.01)
        ax.set_title(og,loc="left")
    fig.tight_layout()
    plt.savefig(fignam, dpi=300, bbox_inches="tight")
    plt.show()


def plot_violin(wd, df):
    # Boxplot
    # f, ax = plt.subplots(3, 4, figsize=(12,8), sharex=True, sharey=True)
    f, ax = plt.subplots(figsize=(6,6))
    month_names = [
                'Measurements','DayCent',
                'DNDC','DNDC-daily'
                ]
    # plot. Set color of marker edge
    flierprops = dict(
                    marker='o', 
                    markerfacecolor='#fc0384', 
                    markersize=7,
                    # linestyle='None',
                    # markeredgecolor='none',
                    alpha=0.3)
    # ax.boxplot(data, flierprops=flierprops)
    r = ax.violinplot(
        df.values,  widths=0.5, showmeans=True, showextrema=True, showmedians=False,
        quantiles=[[0.25, 0.75]]*4,
        bw_method='silverman'
        )
    r['cmeans'].set_color('r')
    r['cquantiles'].set_color('r')
    r['cquantiles'].set_linestyle(':')
    # r['cquantiles'].set_linewidth(3)
    colors = ["#817f82", "#04b0db", '#038f18', '#c40243']
    for c, pc in zip(colors, r['bodies']):
        pc.set_facecolor(c)
    #     pc.set_edgecolor('black')
        pc.set_alpha(0.4)

    ax.set_xticks([i+1 for i in range(4)])
    # ax.set_xticklabels(df_m.keys(), rotation=90)
    ax.set_xticklabels(month_names)
    ax.tick_params(axis='both', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    plt.tight_layout()
    ax.set_ylabel('CH$_4$ emission $(g\;CH_{4}-C\; m^{-2}\cdot d^{-1})$', fontsize=14)
    plt.savefig(os.path.join(wd, 'viloin_plot.png'), dpi=300, bbox_inches="tight")
    plt.show()


def plot_oot(wd, df):

    fig, ax = plt.subplots(figsize=(6,5))
    # colors = cm.tab20(np.linspace(0, 1, len(df.site_name.unique())))
    fmax = df.loc[:, ["ch4_obd", "dc", "dndc", "dndc_d"]].max().max()
    fmin = df.loc[:, ["ch4_obd", "dc", "dndc", "dndc_d"]].min().min()
    colors = ["#04b0db", '#038f18', '#c40243']
    rsqs = []
    for i, c in zip(["dc", "dndc", "dndc_d"], colors):
        x_val = df.loc[:, i].tolist()
        y_val = df.loc[:, "ch4_obd"].tolist()
        correlation_matrix = np.corrcoef(x_val, y_val)
        correlation_xy = correlation_matrix[0,1]
        r_squared = correlation_xy**2
        m, b = np.polyfit(x_val, y_val, 1)
        ax.plot(np.array(x_val), (m*np.array(x_val)) + b, c=c, label='_nolegend_')

        # ax.text(
        #         0.95, 0.05,
        #         f'$y={m:.2f}x{b:.2f}$',
        #         horizontalalignment='right',
        #         # bbox=dict(facecolor='gray', alpha=0.2),
        #         transform=ax.transAxes
        #         )
        ax.scatter(df[i], df.ch4_obd,  alpha=0.4, c=c)
        rsqs.append(r_squared)


    ax.text(
            0.05, 0.9,
            f'$R^2:$ {rsqs[0]:.3f}  |  {rsqs[1]:.3f}  |  {rsqs[2]:.3f}',
            horizontalalignment='left',
            bbox=dict(facecolor='gray', alpha=0.2),
            transform=ax.transAxes
            )





    # rsq_val = round(m1.rsq(df[simnam], df.ch4_obd), 3)
    # rmse_val = round(m1.rmse(df[simnam].values, df.ch4_obd.values), 3)
    # pbias_val = round(m1.pbias(df[simnam].values, df.ch4_obd.values), 3)

    # lgds = []
    # for tn, c in zip(df.site_name.unique(), colors):
    #     sdf = df.loc[df['site_name'] == tn]
    #     ax.scatter(
    #         sdf.somsc_sim, sdf.obd, 
    #         color = c, 
    #         alpha=0.7)
    #     rsq_val = round(rsquared(sdf.obd, sdf.somsc_sim), 3)
    #     rmse_val = round(rmse(sdf.obd, sdf.somsc_sim), 3)
    #     lgds.append(f"{tn} (rsq:{rsq_val}, rmse:{rmse_val})")
    ax.plot([fmin, fmax], [fmin, fmax], 'k--', alpha=0.2)
    ax.set_xlabel("Simulated CH$_4$ emission $(g\;CH_{4}\; m^{-2}\cdot d^{-1})$")
    ax.set_ylabel("Observed CH$_4$ emission $(g\;CH_{4}\; m^{-2}\cdot d^{-1})$")
    # plt.legend(
    #     lgds, 
    #     bbox_to_anchor=(1.05, 1.05), ncols=numcols, fontsize=fsize)
    # fig.tight_layout()
    plt.savefig(os.path.join(wd, "plot_oneToOne.jpg"), dpi=300, bbox_inches="tight")
    plt.show()

def viz(wd):
    m1 = DCmodel(wd)
    obd =  m1.read_inputs()
    obd["date"] = obd.index
    obd["new_idx"] = obd.loc[:, "date"].astype(str) + "-"+ obd.loc[:, "condition"]
    # sim = pd.read_csv(
    #     os.path.join(wd, "ch4_output.csv"), index_col=0, parse_dates=True,)
    simm = pd.read_csv(
        os.path.join(wd, "ch4_multi.out"), sep=r"\s+", comment="#")
    simm["new_idx"] = simm.loc[:, "date"] + "-"+ simm.loc[:, "cont"]
    so_df = simm.merge(obd, how='inner', on='new_idx')
    # so_df = pd.concat([simm.ch4prod, obd.ch4_obd], axis=1).dropna(axis=0)
    os.chdir(wd)
    # plot_one_one(so_df)
    plot_oom(so_df)



