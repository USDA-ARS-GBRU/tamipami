import argparse
import json

def cutoff_arg_validator(value: str) -> dict:
    """
    Validates --cutoff input for tamipami CLI.

    Input must:
    - Be JSON (e.g. '{"3": 0.7, "4": 0.75, ...}')
    - Be a dict with integer keys from 3 to 8 (inclusive)
    - Have values convertible to float

    Args:
        value (str): The --cutoff value as passed via CLI.

    Returns:
        dict[int, float]: Validated and normalized cutoff dictionary.

    Raises:
        argparse.ArgumentTypeError: If input is not valid.
    """
    try:
        data = json.loads(value)
    except Exception:
        raise argparse.ArgumentTypeError(
            "--cutoff must be a valid JSON string, e.g. '{\"3\": 0.7, \"4\": 0.92}'"
        )

    if not isinstance(data, dict):
        raise argparse.ArgumentTypeError("--cutoff must be a JSON object (dictionary)")

    result = {}
    for raw_key, raw_val in data.items():
        # Key as int, in allowed range
        try:
            k = int(raw_key)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Key '{raw_key}' is not an integer. All keys must be integers from 3 to 8."
            )
        if k < 3 or k > 8:
            raise argparse.ArgumentTypeError(
                f"Key '{k}' is outside the valid range (3-8)."
            )
        # Value as float
        try:
            v = float(raw_val)
        except Exception:
            raise argparse.ArgumentTypeError(
                f"Value for key '{k}' must be a number (can be converted to float). Got '{raw_val}'."
            )
        result[k] = v
    return result