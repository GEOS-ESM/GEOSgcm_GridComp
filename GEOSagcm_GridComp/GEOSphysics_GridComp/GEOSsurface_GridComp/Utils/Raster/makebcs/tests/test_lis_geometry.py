#!/usr/bin/env python3

"""Fast regression tests. Does not read the global Discover masks."""

import itertools
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent

sys.path.insert(0, str(SOURCE))

import make_bcs_lis as driver
import make_bcs_lis_tools as tools


def write_raster(path, data, bad_record=False):
    with Path(path).open("wb") as handle:
        nx = data.shape[1]

        for row in data:
            record_size = nx * 4
            header_size = record_size + (4 if bad_record else 0)

            np.array(
                [header_size],
                dtype="<i4",
            ).tofile(handle)

            np.asarray(
                row,
                dtype="<i4",
            ).tofile(handle)

            np.array(
                [record_size],
                dtype="<i4",
            ).tofile(handle)


def write_tile(path, array, im=2, jm=1):
    with Path(path).open("w") as handle:
        handle.write(
            f"{len(array)} {im * jm} 4 2\n"
            "2\n"
            "CF0001x6C\n"
            "1\n"
            "6\n"
            "DE0002xPE0001-LIS\n"
            "2\n"
            "1\n"
        )

        for tile in array:
            handle.write(
                " ".join(
                    str(tile[name])
                    for name in array.dtype.names
                )
                + "\n"
            )


class PlannerTests(unittest.TestCase):

    def setUp(self):
        self.args = driver.parser().parse_args([
            "--resolution",
            "90",
            "--outdir",
            "out",
            "--bin-dir",
            "/bin",
        ])

        self.native = {
            "N_lon": 129600,
            "N_lat": 64800,
            "east_west": 1440,
            "north_south": 720,
            "N_LIS_IDs": 1036800,
        }

    def test_c90_plan(self):
        name, nx, ny, steps = driver.make_plan(
            self.args,
            self.native,
        )

        self.assertEqual(
            name,
            "CF0090x6C_DE1440xPE0720-LIS",
        )

        self.assertEqual(
            (nx, ny),
            (43200, 21600),
        )

        self.assertEqual(len(steps), 5)
        self.assertIn("-l", steps[1][1])
        self.assertIn("1036803", steps[1][1])

        joined = " ".join(
            value
            for _, step in steps
            for value in step
        )

        self.assertNotIn("mkCatchParam", joined)
        self.assertNotIn("Pfafstetter", joined)
        self.assertNotIn("runoff", joined)

    def test_additional_cf_resolutions(self):
        for resolution in (180, 360):
            with self.subTest(resolution=resolution):
                self.args.resolution = resolution

                name = driver.make_plan(
                    self.args,
                    self.native,
                )[0]

                self.assertTrue(
                    name.startswith(
                        f"CF{resolution:04d}x6C_"
                    )
                )

    def test_cube_ocean_plan(self):
        self.args.ocean_grid = "cube"

        name, _, _, steps = driver.make_plan(
            self.args,
            self.native,
        )

        self.assertEqual(
            name,
            "CF0090x6C_CF0090x6C-LIS",
        )

        self.assertEqual(len(steps), 4)

        self.assertEqual(
            steps[-2][1][-2:],
            ["CF0090x6C", "LIS"],
        )

    def test_too_fine_counter_rejected(self):
        self.args.raster_arcsec = 10

        with self.assertRaisesRegex(
            ValueError,
            "counter",
        ):
            driver.make_plan(
                self.args,
                self.native,
            )

    def test_non_nested_LIS_grid_rejected(self):
        self.native["east_west"] = 1436

        with self.assertRaisesRegex(
            ValueError,
            "nest",
        ):
            driver.make_plan(
                self.args,
                self.native,
            )

    def test_nondivisible_ocean_rejected(self):
        self.args.ocean_im = 1439

        with self.assertRaisesRegex(
            ValueError,
            "ocean",
        ):
            driver.make_plan(
                self.args,
                self.native,
            )

    def test_reynolds_reference_not_silently_used(self):
        self.args.reference_mask = Path(
            "GEOS5_10arcsec_mask_freshwater-lakes.nc"
        )

        with self.assertRaisesRegex(
            ValueError,
            "Reynolds",
        ):
            driver.make_plan(
                self.args,
                self.native,
            )

    def test_run_and_dry_run_exclusive(self):
        with self.assertRaises(SystemExit):
            driver.parser().parse_args([
                "--resolution",
                "90",
                "--outdir",
                "out",
                "--run",
                "--dry-run",
            ])

    def test_failure_does_not_report_success(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work = Path(tempdir)
            (work / "logs").mkdir()

            with self.assertRaisesRegex(
                RuntimeError,
                "failed",
            ):
                driver.run_step(
                    "bad",
                    [
                        sys.executable,
                        "-c",
                        "raise SystemExit(17)",
                    ],
                    work,
                    {},
                )

class GeometryTests(unittest.TestCase):

    def fixture(self, work):
        # Four equal-area global quadrants, two native LIS cells.
        # Western cell is land; eastern cell is ocean in both hemispheres.

        tiles = np.zeros(4, dtype=tools.TILE_DTYPE)
        tiles["typ"] = [100, 0, 100, 0]
        tiles["area"] = math.pi
        tiles["lon"] = [-90, 90, -90, 90]
        tiles["lat"] = [-45, -45, 45, 45]
        tiles["i2"] = [1, 2, 1, 2]
        tiles["j2"] = 1

        row = np.array(
            [[1, 1, 2, 2], [3, 3, 4, 4]],
            dtype="i4",
        )

        base = np.zeros(2, dtype=np.dtype(tools.TILE_DTYPE.descr[:8]))
        base["typ"] = [100, 0]
        base["i1"] = [1, 0]
        base["area"] = 2 * math.pi

        base_raster = np.array(
            [[1, 1, 2, 2], [1, 1, 2, 2]],
            dtype="i4",
        )

        write_raster(work / "final.rst", row)
        write_raster(work / "base.rst", base_raster)

        header = {"nx": 4, "ny": 2}

        return header, tiles, base

    def test_raster_record_read(self):
        with tempfile.TemporaryDirectory() as tempdir:
            path = Path(tempdir) / "x.rst"
            expected = np.arange(12).reshape(3, 4)

            write_raster(path, expected)

            actual = np.array(list(tools.raster_rows(path, 4, 3)))
            self.assertTrue(np.array_equal(actual, expected))

    def test_bad_record_rejected(self):
        with tempfile.TemporaryDirectory() as tempdir:
            path = Path(tempdir) / "x.rst"
            write_raster(path, np.zeros((2, 4)), bad_record=True)

            with self.assertRaisesRegex(ValueError, "header"):
                list(tools.raster_rows(path, 4, 2))

    def test_extra_record_rejected(self):
        with tempfile.TemporaryDirectory() as tempdir:
            path = Path(tempdir) / "x.rst"
            write_raster(path, np.zeros((3, 4)))

            with self.assertRaisesRegex(ValueError, "trailing"):
                list(tools.raster_rows(path, 4, 2))

    def test_bounds_and_areas(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work = Path(tempdir)
            header, tiles, base = self.fixture(work)

            bounds, checks = tools.check_raster_and_bounds(
                work / "final.rst",
                work / "base.rst",
                header,
                tiles,
                base,
                2,
                1,
            )

            expected_bounds = [
                [-180, 0, -90, 0],
                [0, 180, -90, 0],
                [-180, 0, 0, 90],
                [0, 180, 0, 90],
            ]

            np.testing.assert_allclose(bounds, expected_bounds)
            self.assertEqual(checks["working_pixels_checked"], 8)

    def test_wrong_LIS_id_rejected(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work = Path(tempdir)
            header, tiles, base = self.fixture(work)
            tiles["i2"][0] = 2

            with self.assertRaisesRegex(ValueError, "longitude assignment"):
                tools.check_raster_and_bounds(
                    work / "final.rst",
                    work / "base.rst",
                    header,
                    tiles,
                    base,
                    2,
                    1,
                )

    def test_wrong_surface_rejected(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work = Path(tempdir)
            header, tiles, base = self.fixture(work)
            tiles["typ"][1] = 19

            with self.assertRaisesRegex(ValueError, "Surface type"):
                tools.check_raster_and_bounds(
                    work / "final.rst",
                    work / "base.rst",
                    header,
                    tiles,
                    base,
                    2,
                    1,
                )

    def test_wrong_area_rejected(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work = Path(tempdir)
            header, tiles, base = self.fixture(work)
            tiles["area"][0] *= 1.01

            with self.assertRaisesRegex(ValueError, "areas"):
                tools.check_raster_and_bounds(
                    work / "final.rst",
                    work / "base.rst",
                    header,
                    tiles,
                    base,
                    2,
                    1,
                )

    def test_ascii_round_trip(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work = Path(tempdir)
            _, tiles, _ = self.fixture(work)

            tiles["i1"] = 1
            tiles["j1"] = 1
            tiles["frac1"] = 1
            tiles["frac2"] = 1

            path = work / "x.til"
            write_tile(path, tiles)

            read_header, read_tiles = tools.read_tiling(path)

            self.assertEqual(read_header["nsource"], 2)
            np.testing.assert_array_equal(read_tiles, tiles)

    def test_mapping_keeps_native_id_and_distinct_weights(self):
        with tempfile.TemporaryDirectory() as tempdir:
            work = Path(tempdir)
            _, tiles, _ = self.fixture(work)

            tiles["i1"] = 1
            tiles["j1"] = [1, 2, 4, 5]
            tiles["frac1"] = 1
            tiles["frac2"] = 0.5

            mapping_path = work / "map.csv"
            tools.export_mapping(mapping_path, tiles, 2, 1, 1)

            rows = mapping_path.read_text().splitlines()

            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[1].split(",")[:4], ["1", "1", "1", "1"])
            self.assertEqual(rows[2].split(",")[:4], ["3", "1", "1", "1"])

    def test_atmospheric_and_LIS_weights(self):
        tiles = np.zeros(6, dtype=tools.TILE_DTYPE)

        tiles["typ"] = 100
        tiles["area"] = 4 * math.pi / 6
        tiles["lon"] = 0
        tiles["lat"] = 0
        tiles["i1"] = 1
        tiles["j1"] = np.arange(1, 7)
        tiles["i2"] = 1
        tiles["j2"] = 1
        tiles["frac1"] = 1
        tiles["frac2"] = 1 / 6

        base = np.zeros(1, dtype=np.dtype(tools.TILE_DTYPE.descr[:8]))
        base["typ"] = 100
        base["i1"] = 1
        base["area"] = 4 * math.pi

        header = {
            "nx": 4,
            "ny": 2,
            "nsource": 2,
            "grids": [("CF0001x6C", 1, 6), ("LIS", 2, 1)],
        }

        result = tools.check_tables(header, tiles, header, base, 2, 1, 1)

        self.assertEqual(result["represented_LIS_cells"], 1)

        tiles["frac1"][0] = 0.8

        with self.assertRaisesRegex(ValueError, "Atmospheric fractions"):
            tools.check_tables(header, tiles, header, base, 2, 1, 1)


@unittest.skipUnless(
    shutil.which("gfortran"),
    "gfortran is needed for compiled geometry tests",
)
class FortranGeometryTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.work = Path(cls.temp.name)

        test_source = cls.work / "test.F90"

        test_source.write_text(
            """program test
  use LISRasterGeometryMod
  implicit none

  integer :: n, i, block(3,3), v

  read(*,*) n

  do i = 1, n
    read(*,*) block
    v = reference_mask_mode(block, 291284)
    write(*,*) v
  end do

  if (lis_id_at_pixel(1, 1, 43200, 21600, 1440, 720) /= 1) stop 1
  if (lis_id_at_pixel(43200, 21600, 43200, 21600, 1440, 720) /= 1036800) stop 2

  if (mask_surface_type(291285, 1036800) /= 100) stop 3
  if (mask_surface_type(291286, 1036800) /= 100) stop 4
  if (mask_surface_type(1036800, 1036800) /= 100) stop 5
  if (mask_surface_type(190000000, 1036800) /= 19) stop 6
  if (mask_surface_type(200000000, 1036800) /= 20) stop 7

end program
"""
        )

        cls.exe = cls.work / "test.x"

        subprocess.run(
            [
                "gfortran", "-O2", "-fcheck=all", "-Wall", "-Wextra",
                str(SOURCE / "LISRasterGeometryMod.F90"),
                str(test_source),
                "-o", str(cls.exe),
            ],
            cwd=cls.work,
            check=True,
            capture_output=True,
            text=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_blocks(self, blocks):
        input_text = (
            str(len(blocks))
            + "\n"
            + "\n".join(" ".join(map(str, block)) for block in blocks)
            + "\n"
        )

        result = subprocess.run(
            [str(self.exe)],
            input=input_text,
            text=True,
            capture_output=True,
            check=True,
        )

        return np.array([int(value) for value in result.stdout.split()])

    def test_exhaustive_2x2_and_random_3x3_vs_legacy_histogram(self):
        states = [0, 1, 2, 190000000, 200000000]

        blocks = [
            list(block) + [0] * 5
            for block in itertools.product(states, repeat=4)
        ]

        rng = np.random.default_rng(1234)
        blocks.extend(rng.choice(states, (4000, 9)).tolist())

        expected = []

        for block in blocks:
            # Reproduce the original temporary IDs, float32 histogram,
            # sorting, and first-maximum behavior. All reference IDs
            # used here are exactly representable as float32 integers.

            values = np.array(block, dtype=np.int64)

            values[values == 190000000] = 291285
            values[values == 200000000] = 291286

            positive = np.unique(values[values > 0])

            if len(positive) == 0:
                expected.append(0)
                continue

            density = np.array([
                (values.astype("f4") == float(value)).sum()
                for value in positive
            ])

            winner = int(positive[density.argmax()])
            winner = {
                291285: 190000000,
                291286: 200000000,
            }.get(winner, winner)

            expected.append(winner)

        np.testing.assert_array_equal(self.run_blocks(blocks), expected)

    def test_surface_not_changed_by_coalescing_LIS_ids(self):
        reference = [1, 1, 1, 2, 2, 190000000, 190000000, 190000000, 190000000]

        self.assertEqual(
            self.run_blocks([reference])[0],
            190000000,
        )

    def test_smallest_id_wins_tie(self):
        blocks = [
            [2, 1, 0, 0, 0, 0, 0, 0, 0],
            [190000000, 1, 0, 0, 0, 0, 0, 0, 0],
            [200000000, 190000000, 0, 0, 0, 0, 0, 0, 0],
        ]

        np.testing.assert_array_equal(
            self.run_blocks(blocks),
            [1, 1, 190000000],
        )


if __name__ == "__main__":
    unittest.main()
