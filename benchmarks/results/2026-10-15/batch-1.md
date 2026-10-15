Every tool built for `native` (Rust `-C target-cpu=native`, C and C++ `-march=native`, Mojo its default, the host), run on Intel(R) Core(TM) i9-7900X CPU @ 3.30GHz.

| workload | threads | dinara-align | dinara-align GPU | WFA2-lib | KSW2 | parasail | Edlib |
| :-- | :-- | --: | --: | --: | --: | --: | --: |
| illumina-affine | one core | **229 ms** | — | 1060 ms | 10.9 s | 21.8 s | — |
| illumina-edit | one core | **211 ms** | — | 448 ms | — | — | 2004 ms |
| illumina-affine | 10 cores | **26 ms** | 34 ms | 112 ms | 1192 ms | 2345 ms | — |
| illumina-edit | 10 cores | **25 ms** | — | 49 ms | — | — | 207 ms |
