from f90nml import Namelist
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.convection.GF_2020.GF_functions.cup_up_moisture import (
    CupUpMoisture,
)


class TranslateCupUpMoisture(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "c1d": {},
            "ccn": {},
            "cd": {},
            "cnvfrc": {},
            "cumulus": {},
            "dby": {},
            "entr_rate": {},
            "gamma_cup": {},
            "hc": {},
            "ierr": {},
            "k22": {},
            "kbcon": {},
            "klcl": {},
            "ktop": {},
            "p_cup": {},
            "po": {},
            "q": {},
            "qe_cup": {},
            "qes_cup": {},
            "rho": {},
            "srftype": {},
            "start_level": {},
            "t_cup": {},
            "up_massdetr": {},
            "up_massentr": {},
            "use_linear_subcl_mf": {},
            "vvel1d": {},
            "vvel2d": {},
            "x_add_buoy": {},
            "xland": {},
            "z_cup": {},
            "zqexec": {},
            "zu": {},
            "zws": {},
            "c0_shal": {},
            "c0_mid": {},
            "c0_deep": {},
            "qrc_crit_lnd": {},
            "qrc_crit_ocn": {},
            "bc_meth": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "clw_all": self.grid.compute_dict(),
            "ierr": self.grid.compute_dict(),
            "psum": self.grid.compute_dict(),
            "psumh": self.grid.compute_dict(),
            "pw": self.grid.compute_dict(),
            "pwav": self.grid.compute_dict(),
            "qc": self.grid.compute_dict(),
            "qrc": self.grid.compute_dict(),
            "tempc": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        cup_up_moisture = CupUpMoisture(
            self.stencil_factory,
            self.grid.quantity_factory,
        )

        # Field inputs

        bc_meth = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(bc_meth.view[:, :, :], inputs["bc_meth"])

        c1d = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(c1d.view[:, :, :], inputs["c1d"])

        ccn = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(ccn.view[:, :, :], inputs["ccn"])

        cd = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(cd.view[:, :, :], inputs["cd"])

        cnvfrc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(cnvfrc.view[:, :, :], inputs["cnvfrc"])

        name = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(name.view[:, :, :], inputs["cumulus"])

        dby = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(dby.view[:, :, :], inputs["dby"])

        entr_rate = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(entr_rate.view[:, :, :], inputs["entr_rate"])

        gamma_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(gamma_cup.view[:, :, :], inputs["gamma_cup"])

        hc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(hc.view[:, :, :], inputs["hc"])

        ierr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(ierr.view[:, :, :], inputs["ierr"])

        k22 = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(k22.view[:, :, :], inputs["k22"])

        kbcon = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(kbcon.view[:, :, :], inputs["kbcon"])

        klcl = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(klcl.view[:, :, :], inputs["klcl"])

        ktop = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(ktop.view[:, :, :], inputs["ktop"])

        p_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(p_cup.view[:, :, :], inputs["p_cup"])

        po = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(po.view[:, :, :], inputs["po"])

        q = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(q.view[:, :, :], inputs["q"])

        qe_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qe_cup.view[:, :, :], inputs["qe_cup"])

        qes_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qes_cup.view[:, :, :], inputs["qes_cup"])

        rho = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(rho.view[:, :, :], inputs["rho"])

        srftype = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(srftype.view[:, :, :], inputs["srftype"])

        start_level = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(start_level.view[:, :, :], inputs["start_level"])

        t_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(t_cup.view[:, :, :], inputs["t_cup"])

        up_massdetr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(up_massdetr.view[:, :, :], inputs["up_massdetr"])

        up_massentr = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(up_massentr.view[:, :, :], inputs["up_massentr"])

        use_linear_subcl_mf = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
            dtype=Int,
        )
        safe_assign_array(
            use_linear_subcl_mf.view[:, :, :], inputs["use_linear_subcl_mf"]
        )

        vvel2d = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(vvel2d.view[:, :, :], inputs["vvel2d"])

        vvel1d = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(vvel1d.view[:, :, :], inputs["vvel1d"])

        x_add_buoy = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(x_add_buoy.view[:, :, :], inputs["x_add_buoy"])

        xland = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(xland.view[:, :, :], inputs["xland"])

        z_cup = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(z_cup.view[:, :, :], inputs["z_cup"])

        zqexec = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(zqexec.view[:, :, :], inputs["zqexec"])

        zu = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(zu.view[:, :, :], inputs["zu"])

        zws = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(zws.view[:, :, :], inputs["zws"])

        c0_shal = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(c0_shal.view[:, :, :], inputs["c0_shal"])

        c0_mid = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(c0_mid.view[:, :, :], inputs["c0_mid"])

        c0_deep = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(c0_deep.view[:, :, :], inputs["c0_deep"])

        qrc_crit_lnd = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qrc_crit_lnd.view[:, :, :], inputs["qrc_crit_lnd"])

        qrc_crit_ocn = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )
        safe_assign_array(qrc_crit_ocn.view[:, :, :], inputs["qrc_crit_ocn"])

        clw_all = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        psum = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        psumh = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        pw = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        pw = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        pwav = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        qc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        qrc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        tempc = QuantityFactory.zeros(
            self.quantity_factory,
            dims=[X_DIM, Y_DIM, Z_DIM],
            units="n/a",
        )

        print("running stencil")
        cup_up_moisture(
            # In
            c1d=c1d,
            ccn=ccn,
            cd=cd,
            cnvfrc=cnvfrc,
            name=name,
            dby=dby,
            entr_rate=entr_rate,
            gamma_cup=gamma_cup,
            hc=hc,
            k22=k22,
            kbcon=kbcon,
            klcl=klcl,
            ktop=ktop,
            p_cup=p_cup,
            po=po,
            q=q,
            qe_cup=qe_cup,
            qes_cup=qes_cup,
            rho=rho,
            srftype=srftype,
            start_level=start_level,
            t_cup=t_cup,
            up_massdetr=up_massdetr,
            up_massentr=up_massentr,
            use_linear_subcl_mf=use_linear_subcl_mf,
            vvel1d=vvel1d,
            vvel2d=vvel2d,
            x_add_buoy=x_add_buoy,
            xland=xland,
            z_cup=z_cup,
            zqexec=zqexec,
            zu=zu,
            zws=zws,
            c0_shal=c0_shal,
            c0_mid=c0_mid,
            c0_deep=c0_deep,
            qrc_crit_lnd=qrc_crit_lnd,
            qrc_crit_ocn=qrc_crit_ocn,
            bc_meth=bc_meth,
            # Out
            clw_all=clw_all,
            ierr=ierr,
            psum=psum,
            psumh=psumh,
            pw=pw,
            pwav=pwav,
            qc=qc,
            qrc=qrc,
            tempc=tempc,
        )

        return {
            "clw_all": clw_all.view[:],
            "ierr": ierr.view[:],
            "psum": psum.view[:],
            "psumh": psumh.view[:],
            "pw": pw.view[:],
            "pwav": pwav.view[:],
            "qc": qc.view[:],
            "qrc": qrc.view[:],
            "tempc": tempc.view[:],
        }
