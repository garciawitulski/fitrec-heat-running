"""Recreate the paper's random frequencies. Inputs are counts, never runner IDs.

The reference used numpy 2.2.6, PCG64 via default_rng. R sorts the runner roster
and applies the integer frequencies, fits all models and computes the intervals.
The binary file is a local temporary exchange, with little-endian int32 counts.
"""
import argparse

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path")
    parser.add_argument("design", choices=["joint", "lee"])
    parser.add_argument("sizes", nargs="+", type=int)
    args = parser.parse_args()
    with open(args.path, "wb") as output:
        if args.design == "joint":
            (size,) = args.sizes
            for replicate in range(999):
                rng = np.random.default_rng([20261002, replicate])
                draw = rng.integers(0, size, size)
                np.bincount(draw, minlength=size).astype("<i4").tofile(output)
        else:
            rng = np.random.default_rng(20261104)
            # Outcome order: whole-session heart rate, speed. Within each outcome:
            # point-estimate trim, lower-limit trim. The stream advances across all.
            for size in args.sizes:
                for _ in range(2 * 199):
                    draw = rng.integers(0, size, size)
                    np.bincount(draw, minlength=size).astype("<i4").tofile(output)


if __name__ == "__main__":
    main()
