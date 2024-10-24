from . import main
import math
import statistics

class TestCheckN:

    # Calculate fraction of SD for a list of positive integers
    def test_fraction_of_sd_positive_integers(self):
        vect = [4, 9, 16, 25]
        result = main.check_N(vect)
        expected = math.sqrt(statistics.mean(vect)) / statistics.stdev(vect)
        assert math.isclose(result, expected, rel_tol=1e-9)
    # Handle an empty list without errors
    def test_handle_empty_list(self):
        vect = []
        with pytest.raises(statistics.StatisticsError):
            main.check_N(vect)