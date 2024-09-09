import ndsl.stencils.testing.conftest
import savepoint
from ndsl.stencils.testing.conftest import *  # noqa: F403,F401

# Point to an __init__.py where all the TestX are improted
ndsl.stencils.testing.conftest.translate = savepoint  # type: ignore
