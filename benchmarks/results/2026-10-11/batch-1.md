Every tool built for `native` (Rust `-C target-cpu=native`, C and C++ `-march=native`, Mojo its default, the host), run on Intel(R) Core(TM) i9-7900X CPU @ 3.30GHz.

| workload | threads | dinara-align | dinara-align GPU | WFA2-lib | KSW2 | parasail | Edlib |
| :-- | :-- | --: | --: | --: | --: | --: | --: |
| illumina-affine | one core | **231 ms** | — | 1063 ms | 10.9 s | 21.7 s | — |
| illumina-edit | one core | **212 ms** | — | 450 ms | — | — | 1971 ms |
| illumina-affine | 10 cores | **30 ms** | 49 ms | 113 ms | 1198 ms | 2350 ms | — |
| illumina-edit | 10 cores | **27 ms** | — | 56 ms | — | — | 213 ms |
