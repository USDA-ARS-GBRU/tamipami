import pytest

from tamipami import app


class TestReadCountCheck:
    """Checks that read_count_check flags runs where too few reads hit a target.

    read_count_check calls st.warning when the fraction of reads containing a
    guide target drops below minfrac, so the tests patch app.st.warning and
    look at whether it was called.
    """

    def test_low_target_fraction_fires_warning(self, mocker):
        warn = mocker.patch.object(app.st, "warning")
        run_summ = {
            "exp": {"tot": 1000, "targets": 50},
            "cont": {"tot": 1000, "targets": 40},
        }
        app.read_count_check(run_summ, minfrac=0.5)
        warn.assert_called_once()

    def test_bad_control_only_fires_warning(self, mocker):
        # Experimental fraction is fine (0.9) but the control fraction is bad
        # (0.05). Both sides must be checked, so the warning should still fire.
        warn = mocker.patch.object(app.st, "warning")
        run_summ = {
            "exp": {"tot": 1000, "targets": 900},
            "cont": {"tot": 1000, "targets": 50},
        }
        app.read_count_check(run_summ, minfrac=0.5)
        warn.assert_called_once()

    def test_healthy_run_is_quiet(self, mocker):
        warn = mocker.patch.object(app.st, "warning")
        run_summ = {
            "exp": {"tot": 1000, "targets": 800},
            "cont": {"tot": 1000, "targets": 750},
        }
        app.read_count_check(run_summ, minfrac=0.5)
        warn.assert_not_called()
