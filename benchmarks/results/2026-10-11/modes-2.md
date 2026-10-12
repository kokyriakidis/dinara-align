Every tool built for `native` (Rust `-C target-cpu=native`, C and C++ `-march=native`, Mojo its default, the host), run on Intel(R) Core(TM) i9-7900X CPU @ 3.30GHz.

| workload | dinara-align | Edlib | WFA2-lib | KSW2 | parasail | SSW |
| :-- | --: | --: | --: | --: | --: | --: |
| infix-edit | **130 µs** | 268 µs | 1 ms | — | — | — |
| infix-edit-long | **2.01 ms** | 13.6 ms | 54 ms | — | — | — |
| prefix-edit | **42.5 µs** | 115 µs | 47 µs | — | — | — |
| infix-affine | **3.92 ms** | — | 3.97 ms | — | — | — |
| extension | **726 µs** | — | — | 2.06 ms | — | — |
| extension-bonus | **166 µs** | — | — | 1.24 ms | — | — |
| two-piece | **8.25 ms** | — | 15.9 ms | 46.9 ms | — | — |
| table-global | **196 µs** | — | — | — | 2.94 ms | — |
| table-local | **1.91 ms** | — | — | — | 5.15 ms | 3.28 ms |
| local-short | **32 µs** | — | — | — | 68.2 µs | 51.2 µs |
| local-window | **1.42 ms** | — | — | — | 4.75 ms | 3.15 ms |
| local-long | **17.6 ms** | — | — | — | 293 ms | 84.5 ms |
| overlap-reads | **737 µs** | — | — | — | 15.3 ms | — |
| infix-reads | **687 µs** | — | — | — | 15 ms | — |
