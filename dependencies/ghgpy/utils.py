import numpy as np
import datetime
import os
from ghgpy.models import DCmodel
import pandas as pd

class ObjFns:
    def __init__(self) -> None:
        pass
        
    # def obj_fns(obj_fn, sims, obds):
    #     return obj_fn(sims, obds)
    @staticmethod
    def nse(sims, obds):
        """Nash-Sutcliffe Efficiency (NSE) as per `Nash and Sutcliffe, 1970
        <https://doi.org/10.1016/0022-1694(70)90255-6>`_.

        .. math::
            E_{\\text{NSE}} = 1 - \\frac{\\sum_{i=1}^{N}[e_{i}-s_{i}]^2}
            {\\sum_{i=1}^{N}[e_{i}-\\mu(e)]^2}

        where *N* is the length of the *sims* and *obds*
        periods, *e* is the *obds* series, *s* is (one of) the
        *sims* series, and *μ* is the arithmetic mean.

        """
        nse_ = 1 - (
                np.sum((obds - sims) ** 2, axis=0, dtype=np.float64)
                / np.sum((obds - np.mean(obds)) ** 2, dtype=np.float64)
        )
        return round(nse_, 4)

    @staticmethod
    def rmse(sims, obds):
        """Root Mean Square Error (RMSE).

        .. math::
            E_{\\text{RMSE}} = \\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}[e_i-s_i]^2}

        where *N* is the length of the *sims* and *obds*
        periods, *e* is the *obds* series, *s* is (one of) the
        *sims* series.

        """
        rmse_ = np.sqrt(np.mean((obds - sims) ** 2,
                                axis=0, dtype=np.float64))

        return round(rmse_, 4)

    @staticmethod
    def pbias(sims, obds):
        """Percent Bias (PBias).

        .. math::
            E_{\\text{PBias}} = 100 × \\frac{\\sum_{i=1}^{N}(e_{i}-s_{i})}{\\sum_{i=1}^{N}e_{i}}

        where *N* is the length of the *sims* and *obds*
        periods, *e* is the *obds* series, and *s* is (one of)
        the *sims* series.

        """
        pbias_ = (100 * np.sum(obds - sims, axis=0, dtype=np.float64)
                / np.sum(obds))

        return round(pbias_, 4)

    @staticmethod
    def rsq(sims, obds):
        ## R-squared
        rsq_ = (
            (
                (sum((obds - obds.mean())*(sims-sims.mean())))**2
            ) 
            /
            (
                (sum((obds - obds.mean())**2)* (sum((sims-sims.mean())**2))
            ))
        )
        return round(rsq_, 4)
    
class PostPr(object):
    # get sim CH4 emission and obd
    def __init__(self, wd):
        self.wd = wd
        os.chdir(self.wd)

    def get_ch4_so_df(self, outfnam=None):
        if outfnam is None:
            outfnam = "ch4_multi.out"
        # m1 = DCmodel(self.wd)
        # obd =  m1.read_inputs()
        # obd["date"] = obd.index
        # obd["new_idx"] = obd.loc[:, "date"].astype(str) + "-"+ obd.loc[:, "condition"]
        # sim = pd.read_csv(
        #     os.path.join(wd, "ch4_output.csv"), index_col=0, parse_dates=True,)
        simm = pd.read_csv(
            os.path.join(self.wd, outfnam), sep=r"\s+", comment="#")
        simm["new_idx"] = simm.loc[:, "date"] + "-"+ simm.loc[:, "cont"]
        # so_df = simm.merge(obd, how='inner', on='new_idx')
        return simm

        
if __name__ == "__main__":
    wd = "D:\\Projects\\Models\\swatp-ghg\\models"
    df = PostPr(wd).get_ch4_so_df()
    print(df.columns)
    rrr = ObjFns.rsq(df.ch4e_tot.values, df.ch4_obd.values)
    print(rrr)
