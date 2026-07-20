try:
    import matplotlib as mpl

    # chose the matplotlib backend that can only generate pngs
    # but not interact with the os (better for testing)
    mpl.use("Agg")
except ImportError:
    pass

import pytest


def pytest_addoption(parser):
    """Add a custom CLI option run all tests including those marked as remote."""
    parser.addoption(
        "--run-remote",
        action="store_true",
        default=False,
        help="Run all tests regardless of markers",
    )


# function executed right after test items collected but before test run
def pytest_collection_modifyitems(config, items):
    """Add option to run remote_data tests or not."""
    if not config.getoption("-m"):
        skip_me = pytest.mark.skip(reason="use `-m remote_data` to run this test")
        for item in items:
            if "remote_data" in item.keywords and not config.getoption("--run-remote"):
                item.add_marker(skip_me)
