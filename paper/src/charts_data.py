"""Figure data, from benchmarks/README.md's unit-cost run over QuickEd's simulated sets (2026-10-11, the
Skylake-X, dinara-align at `dc980d7`, a 30 s budget). Times in ms per alignment.

A value is [ms] or [ms, 1] where the tool finished only some pairs; None where it finished none."""

DIVERGENCE_X = [1, 5, 10, 20]
DIVERGENCE = {
    "dinara-align": [[1.5], [5.02], [7.21], [16.0]],
    "A*PA2-full": [[7.53], [8.45], [13.9], [36.7]],
    "A*PA2-simple": [[3.26], [14.9], [27.6], [56.8]],
    "A*PA": [[7.76], [14.4], [48.3], None],
    "Edlib": [[57.6], [146.0], [247.0], [446.0]],
    "BiWFA": [[3.71], [63.3], [218.0], [695.0]],
    "WFA": [[4.05], [55.0], [182.0], [551.0]],
    "QuickEd": [[12.5], [24.0], [39.0], [62.8]]}
LENGTH_X = [10, 100, 500, 1000]
LENGTH = {
    5: {
        "dinara-align": [[0.269], [5.02], [31.0], [64.1]],
        "A*PA2-full": [[0.797], [8.45], [45.0], [118.0]],
        "A*PA2-simple": [[0.288], [14.9], [214.0], [818.0]],
        "A*PA": [[1.08], [14.4], [83.1], [183.0]],
        "Edlib": [[3.79], [146.0], [2280.0], [8670.0, 1]],
        "BiWFA": [[0.642], [63.3], [1560.0], [6350.0]],
        "WFA": [[0.508], [55.0], [1480.0], None],
        "QuickEd": [[1.25], [24.0], [401.0], [1440.0]]},
    20: {
         "dinara-align": [[0.519], [16.0], [254.0], [901.0]],
         "A*PA2-full": [[1.44], [36.7], [545.0], [2020.0]],
         "A*PA2-simple": [[1.02], [56.8], [901.0], [3520.0]],
         "A*PA": [[205.0], None, None, None],
         "Edlib": [[5.93], [446.0], [8520.0, 1], None],
         "BiWFA": [[7.33], [695.0], [18000.0, 1], None],
         "WFA": [[5.32], [551.0], None, None],
         "QuickEd": [[1.55], [62.8], [1240.0], [4800.0]]},
}
LENGTH20 = LENGTH[20]


def js():
    """The data as the chart script's constants."""
    import json
    return (f"const DIVERGENCE = {{ x: {json.dumps(DIVERGENCE_X)}, unit: \"%\", series: {json.dumps(DIVERGENCE)} }};\n"
            f"const LENGTH_X = {json.dumps(LENGTH_X)};\n"
            f"const LENGTH = {json.dumps(LENGTH)};\n")
