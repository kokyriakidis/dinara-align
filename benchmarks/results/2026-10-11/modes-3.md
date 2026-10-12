Every tool built for `native` (Rust `-C target-cpu=native`, C and C++ `-march=native`, Mojo its default, the host), run on Intel(R) Core(TM) i9-7900X CPU @ 3.30GHz.

| workload | dinara-align | Edlib | WFA2-lib | KSW2 | parasail | SSW |
| :-- | --: | --: | --: | --: | --: | --: |
| infix-edit | **130 µs** | 268 µs | 1 ms | — | — | — |
| infix-edit-long | **2.01 ms** | 13.3 ms | 53.9 ms | — | — | — |
| prefix-edit | **42.4 µs** | 114 µs | 47 µs | — | — | — |
| infix-affine | **3.92 ms** | — | 3.95 ms | — | — | — |
| extension | **729 µs** | — | — | 2.06 ms | — | — |
| extension-bonus | **166 µs** | — | — | 1.24 ms | — | — |
| two-piece | **8.36 ms** | — | 15.8 ms | 46.9 ms | — | — |
| table-global | **196 µs** | — | — | — | 2.93 ms | — |
| table-local | **1.88 ms** | — | — | — | 5.16 ms | 3.27 ms |
| local-short | **32 µs** | — | — | — | 67.2 µs | 51.3 µs |
| local-window | **1.36 ms** | — | — | — | 4.75 ms | 3.17 ms |
| local-long | **17.4 ms** | — | — | — | 293 ms | 84.8 ms |
| overlap-reads | **702 µs** | — | — | — | 15.3 ms | — |
| infix-reads | **687 µs** | — | — | — | 15.1 ms | — |
