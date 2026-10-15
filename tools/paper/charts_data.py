"""Figure data, from benchmarks/README.md's unit-cost run over QuickEd's simulated sets (2026-10-11, the
Skylake-X, dinara-align at `036e12b`, a 30 s budget). Times in ms per alignment.

A value is [ms] or [ms, 1] where the tool finished only some pairs; None where it finished none."""

DIVERGENCE_X = [1, 5, 10, 20]
DIVERGENCE = {
    "dinara-align": [[1.11], [4.84], [7.04], [15.4]],
    "A*PA2-full": [[7.53], [8.35], [13.9], [36.8]],
    "A*PA2-simple": [[3.26], [14.9], [27.7], [56.8]],
    "A*PA": [[7.76], [14.1], [48.4], None],
    "Edlib": [[58.1], [147.0], [247.0], [445.0]],
    "BiWFA": [[3.72], [63.3], [218.0], [709.0]],
    "WFA": [[4.08], [54.9], [182.0], [553.0]],
    "QuickEd": [[12.5], [23.9], [39.0], [62.6]]}
LENGTH_X = [10, 100, 500, 1000]
LENGTH = {
    5: {
        "dinara-align": [[0.187], [4.84], [29.7], [61.6]],
        "A*PA2-full": [[0.795], [8.35], [44.7], [116.0]],
        "A*PA2-simple": [[0.288], [14.9], [213.0], [819.0]],
        "A*PA": [[1.08], [14.1], [81.7], [183.0]],
        "Edlib": [[3.81], [147.0], [2290.0], [8830.0, 1]],
        "BiWFA": [[0.637], [63.3], [1550.0], [6350.0]],
        "WFA": [[0.513], [54.9], [1480.0], None],
        "QuickEd": [[1.26], [23.9], [394.0], [1440.0]]},
    20: {
         "dinara-align": [[0.439], [15.4], [251.0], [886.0]],
         "A*PA2-full": [[1.44], [36.8], [544.0], [2020.0]],
         "A*PA2-simple": [[1.01], [56.8], [898.0], [3510.0]],
         "A*PA": [[204.0], None, None, None],
         "Edlib": [[5.96], [445.0], [8530.0, 1], None],
         "BiWFA": [[7.31], [709.0], [18000.0, 1], None],
         "WFA": [[5.34], [553.0], None, None],
         "QuickEd": [[1.54], [62.6], [1240.0], [4800.0]]},
}
LENGTH20 = LENGTH[20]


def js():
    """The data as the chart script's constants."""
    import json
    return (f"const DIVERGENCE = {{ x: {json.dumps(DIVERGENCE_X)}, unit: \"%\", series: {json.dumps(DIVERGENCE)} }};\n"
            f"const LENGTH_X = {json.dumps(LENGTH_X)};\n"
            f"const LENGTH = {json.dumps(LENGTH)};\n")
