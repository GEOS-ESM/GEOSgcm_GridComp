from ndsl import Namelist, StencilFactory, Quantity
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.GFDL_1M.stencils import calculate_derived_states, find_klcl, dumb_stencil
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM


class Translatederived_states(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # grid.compute_dict is workaround to remove grid halo, which is hardcoded to 3
        self.in_vars["data_vars"] = {
            "PLE": grid.compute_dict(),
            "ZLE": grid.compute_dict(),
            "T": grid.compute_dict(),
            "U": grid.compute_dict(),
            "V": grid.compute_dict(),
            "Q": grid.compute_dict(),
        }

        self.out_vars = {
            "PLEmb": grid.compute_dict(),
            "PLmb": grid.compute_dict(),
            "ZLE0": grid.compute_dict(),
            "ZL0": grid.compute_dict(),
            "DZET": grid.compute_dict(),
            "DP": grid.compute_dict(),
            "MASS": grid.compute_dict(),
            "iMASS": grid.compute_dict(),
            "U0": grid.compute_dict(),
            "V0": grid.compute_dict(),
            "QST3": grid.compute_dict(),
            "DQST3": grid.compute_dict(),
            "KLCL": grid.compute_dict(),
            "LTS": grid.compute_dict(),
            "EIS": grid.compute_dict(),
            "EXTRADATA": grid.compute_dict(),
        }

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize extra quantities
        self.temporaries = Temporaries.make(self.quantity_factory)
        self.masks = Masks.make(self.quantity_factory)

        # Construct stencils
        self.dumb_stencil = stencil_factory.from_dims_halo(
            func=dumb_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.calculate_derived_states = stencil_factory.from_dims_halo(
            func=calculate_derived_states,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_klcl = stencil_factory.from_dims_halo(
            func=find_klcl,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

    def make_ijk_quantity(self, data, interface: bool = False) -> Quantity:
        if interface == True:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity
        else:
            quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], "n/a")
            quantity.view[:, :, :] = quantity.np.asarray(data[:, :, :])
            return quantity

    def compute(self, inputs):
        p_interface = self.make_ijk_quantity(inputs.pop("PLE"), True)
        geopotential_height_interface = self.make_ijk_quantity(inputs.pop("ZLE"), interface=True)
        t = self.make_ijk_quantity(inputs.pop("T"))
        u = self.make_ijk_quantity(inputs.pop("U"))
        v = self.make_ijk_quantity(inputs.pop("V"))
        vapor = self.make_ijk_quantity(inputs.pop("Q"))

        # self.dumb_stencil(in_data=u, out_data=self.temporaries.u_unmodified)

        self.calculate_derived_states(
            p_interface=p_interface,
            p_interface_mb=self.temporaries.p_interface_mb,
            p_mb=self.temporaries.p_mb,
            geopotential_height_interface=geopotential_height_interface,
            edge_height_above_surface=self.temporaries.edge_height_above_surface,
            layer_height_above_surface=self.temporaries.layer_height_above_surface,
            layer_thickness=self.temporaries.layer_thickness,
            dp=self.temporaries.dp,
            mass=self.temporaries.mass,
            t=t,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            qsat=self.temporaries.qsat,
            dqsat=self.temporaries.dqsat,
            u=u,
            u_unmodified=self.temporaries.u_unmodified,
            v=v,
            v_unmodified=self.temporaries.v_unmodified,
        )

        self.find_klcl(
            t=t,
            p_mb=self.temporaries.p_mb,
            vapor=vapor,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            found_level=self.masks.boolean_2d_mask,
            k_lcl=self.temporaries.k_lcl,
            test_field=self.temporaries.test_field,
        )

        return {
            "PLEmb": self.temporaries.p_interface_mb.field[:],
            "PLmb": self.temporaries.p_mb.field[:],
            "ZLE0": self.temporaries.edge_height_above_surface.field[:],
            "ZL0": self.temporaries.layer_height_above_surface.field[:],
            "DZET": self.temporaries.layer_thickness.field[:],
            "DQST3": self.temporaries.dqsat.field[:],
            "DP": self.temporaries.dp.field[:],
            "MASS": self.temporaries.mass.field[:],
            "iMASS": 1 / self.temporaries.mass.field[:],
            "U0": self.temporaries.u_unmodified.field[:],
            "V0": self.temporaries.v_unmodified.field[:],
            "QST3": self.temporaries.qsat.field[:],
            "DQST3": self.temporaries.dqsat.field[:],
            "KLCL": self.temporaries.k_lcl.field[:],
            "LTS": self.temporaries.lower_tropospheric_stability.field[:],
            "EIS": self.temporaries.estimated_inversion_strength.field[:],
            "EXTRADATA": self.temporaries.test_field.field[:],
        }
