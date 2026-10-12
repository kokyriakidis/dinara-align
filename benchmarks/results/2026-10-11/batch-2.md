Every tool built for `native` (Rust `-C target-cpu=native`, C and C++ `-march=native`, Mojo its default, the host), run on Intel(R) Core(TM) i9-7900X CPU @ 3.30GHz.

| workload | threads | dinara-align | dinara-align GPU | WFA2-lib | KSW2 | parasail | Edlib |
| :-- | :-- | --: | --: | --: | --: | --: | --: |
| illumina-affine | one core | **230 ms** | — | 1070 ms | 10.9 s | 21.7 s | — |
| illumina-edit | one core | **212 ms** | — | 450 ms | — | — | 2037 ms |
| illumina-affine | 10 cores | **27 ms** | 35 ms | 115 ms | 1207 ms | 2340 ms | — |
| illumina-edit | 10 cores | **29 ms** | — | 47 ms | — | — | 209 ms |
