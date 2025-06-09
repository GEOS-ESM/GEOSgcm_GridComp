from ndsl import Namelist, StencilFactory, Quantity
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.stencils import calculate_derived_states, find_klcl, vertical_interpolation, find_eis
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float


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
        }

        # Initalize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        # Initalize extra quantities
        self.outputs = Outputs.make(self.quantity_factory)
        self.temporaries = Temporaries.make(self.quantity_factory)
        self.masks = Masks.make(self.quantity_factory)

        # Construct stencils
        self.calculate_derived_states = stencil_factory.from_dims_halo(
            func=calculate_derived_states,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_klcl = stencil_factory.from_dims_halo(
            func=find_klcl,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.vertical_interpolation = stencil_factory.from_dims_halo(
            func=vertical_interpolation,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
        )

        self.find_eis = stencil_factory.from_dims_halo(
            func=find_eis,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
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

        self.calculate_derived_states(
            p_interface=p_interface,
            p_interface_mb=self.temporaries.p_interface_mb,
            p_mb=self.temporaries.p_mb,
            geopotential_height_interface=geopotential_height_interface,
            edge_height_above_surface=self.temporaries.edge_height_above_surface,
            layer_height_above_surface=self.temporaries.layer_height_above_surface,
            layer_thickness=self.temporaries.layer_thickness,
            layer_thinkness_negative=self.temporaries.layer_thickness_negative,
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
            temporary_3d=self.temporaries.temporary_3d,
            th=self.temporaries.th,
        )

        self.find_klcl(
            t=t,
            p_mb=self.temporaries.p_mb,
            vapor=vapor,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            found_level=self.masks.boolean_2d_mask,
            k_lcl=self.temporaries.k_lcl,
        )

        self.vertical_interpolation(
            field=self.temporaries.th,
            interpolated_field=self.temporaries.th700,
            p_interface_mb=self.temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=self.temporaries.temporary_2d_1,
            pt=self.temporaries.temporary_2d_2,
            boolean_2d_mask=self.masks.boolean_2d_mask,
        )

        self.vertical_interpolation(
            field=t,
            interpolated_field=self.temporaries.t700,
            p_interface_mb=self.temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=self.temporaries.temporary_2d_1,
            pt=self.temporaries.temporary_2d_2,
            boolean_2d_mask=self.masks.boolean_2d_mask,
        )
        self.vertical_interpolation(
            field=self.temporaries.layer_height_above_surface,
            interpolated_field=self.temporaries.z700,
            p_interface_mb=self.temporaries.p_interface_mb,
            target_pressure=Float(70000.0),
            pb=self.temporaries.temporary_2d_1,
            pt=self.temporaries.temporary_2d_2,
            boolean_2d_mask=self.masks.boolean_2d_mask,
        )

        self.find_eis(
            t=t,
            th=self.temporaries.th,
            layer_height_above_surface=self.temporaries.layer_height_above_surface,
            t700=self.temporaries.t700,
            th700=self.temporaries.th700,
            z700=self.temporaries.z700,
            k_lcl=self.temporaries.k_lcl,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
            lower_tropospheric_stability=self.outputs.lower_tropospheric_stability,
            estimated_inversion_strength=self.outputs.estimated_inversion_strength,
        )

        return {
            "PLEmb": self.temporaries.p_interface_mb.field,
            "PLmb": self.temporaries.p_mb.field,
            "ZLE0": self.temporaries.edge_height_above_surface.field,
            "ZL0": self.temporaries.layer_height_above_surface.field,
            "DZET": self.temporaries.layer_thickness.field,
            "DQST3": self.temporaries.dqsat.field,
            "DP": self.temporaries.dp.field,
            "MASS": self.temporaries.mass.field,
            "iMASS": 1 / self.temporaries.mass.field,
            "U0": self.temporaries.u_unmodified.field,
            "V0": self.temporaries.v_unmodified.field,
            "QST3": self.temporaries.qsat.field,
            "DQST3": self.temporaries.dqsat.field,
            "KLCL": self.temporaries.k_lcl.field + 1,  # add 1 b/c python indexing starts at 0
            "LTS": self.outputs.lower_tropospheric_stability.field,
            "EIS": self.outputs.estimated_inversion_strength.field,
        }
