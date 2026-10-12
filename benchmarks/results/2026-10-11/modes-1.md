Every tool built for `native` (Rust `-C target-cpu=native`, C and C++ `-march=native`, Mojo its default, the host), run on Intel(R) Core(TM) i9-7900X CPU @ 3.30GHz.

| workload | dinara-align | Edlib | WFA2-lib | KSW2 | parasail | SSW |
| :-- | --: | --: | --: | --: | --: | --: |
| infix-edit | **130 µs** | 272 µs | 1.01 ms | — | — | — |
| infix-edit-long | **2.02 ms** | 13.5 ms | 54 ms | — | — | — |
| prefix-edit | **42.4 µs** | 114 µs | 47.2 µs | — | — | — |
| infix-affine | **3.92 ms** | — | 3.96 ms | — | — | — |
| extension | **726 µs** | — | — | 2.06 ms | — | — |
| extension-bonus | **166 µs** | — | — | 1.23 ms | — | — |
| two-piece | **8.25 ms** | — | 15.8 ms | 46.9 ms | — | — |
| table-global | **196 µs** | — | — | — | 2.97 ms | — |
| table-local | **1.93 ms** | — | — | — | 5.16 ms | 3.28 ms |
| local-short | **32.1 µs** | — | — | — | 66.9 µs | 51.1 µs |
| local-window | **1.36 ms** | — | — | — | 4.72 ms | 3.15 ms |
| local-long | **17.3 ms** | — | — | — | 294 ms | 84.7 ms |
| overlap-reads | **703 µs** | — | — | — | 15.3 ms | — |
| infix-reads | **691 µs** | — | — | — | 15.1 ms | — |
