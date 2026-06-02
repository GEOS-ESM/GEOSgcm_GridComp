import os

import translate_tests


os.environ["NDSL_LITERAL_PRECISION"] = "32"

import ndsl.stencils.testing.conftest  # noqa: E402
from ndsl.stencils.testing.conftest import *  # noqa: F403,F401, E402


# Point to an __init__.py where all the TestX are improted
ndsl.stencils.testing.conftest.translate = translate_tests  # type: ignore
