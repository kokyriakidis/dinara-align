# Idiomatic Mojo

The main `mojo-syntax` skill makes code compile; these rules make it pass
review. Before submitting or reviewing Mojo, scan the diff for every WRONG
pattern below. The first five sections are the ones reviewers flag most
often in agent-written code.

## 1. Use existing helpers; never open-code them

Search `std.math`, `std.gpu`, `nn/`, `linalg/`, and the TileTensor API
before writing index math, reductions, reshapes, or activations. The
existing code a few lines away usually already uses the helper.

```mojo
# WRONG                                  # CORRECT
(M + BM - 1) // BM                       ceildiv(M, BM)          # std.math
(x + a - 1) // a * a                     align_up(x, a)          # std.math
var s = b // H; var h = b % H            var s, h = divmod(b, H)
1 / (1 + exp(-x))                        sigmoid(x)              # nn.activations
manual shuffle-reduce loop               warp.sum(val)           # std.gpu warp
hard-coded 32 lanes                      WARP_SIZE
scalar tap loop + cast per element       (v.load[w](i) * wts).reduce_add(), cast after reduce
manual ", " string-building loop         ", ".join(items)
```

`ceildiv` works on `Int` and on any `SIMD` scalar (`UInt32`, ...).
Prefer `vectorize` over hand-written scalar loops, and call LLVM/NVVM
intrinsics rather than inline assembly.

Vectorized pointer loads and stores default to the scalar alignment, which
can split them into scalar accesses. Pass the real alignment; the address
must actually satisfy it, or the access faults:

```mojo
# WRONG
var v = ptr.unsafe_load[width=8](off)
# CORRECT
comptime align = align_of[SIMD[dtype, 8]]()
var v = ptr.unsafe_load[width=8, alignment=align](off)
```

## 2. Don't cast values that already have the right type

`thread_idx`, `block_idx`, `block_dim`, `grid_dim` fields return `Int`.
Don't wrap them in `Int(...)`, and don't copy that pattern from older
code. Keep a value in its natural type; a downcast followed by an upcast
is a review flag.

```mojo
# WRONG                                           # CORRECT
var token = Int(block_idx.x)                      var token = block_idx.x
var tile_m = UInt32(block_idx.y)  ...  Int(tile_m)   var tile_m = block_idx.y
```

Cast only at the boundary where the other type is actually required
(for example, comparing against a `UInt32` offset).

## 3. Let parameters be inferred

Omit parameter lists that the arguments already determine. Don't add a
`rank`/`dtype`/`layout` parameter that the argument type carries.

```mojo
# WRONG
Self.run_at_tile[c_layout, a_layout, b_layout, c_engine, a_engine, b_engine](c, a, b, m, n)
epilogue[dtype, width, alignment=alignment](idx, val)
def f[rank: Int, axis: Int](shape: IndexList[rank]): ...
x: TileTensor[DType.float32, layout, origin, Engine=engine]

# CORRECT
Self.run_at_tile(c, a, b, m, n)
epilogue[alignment=alignment](idx, val)
def f[axis: Int](shape: IndexList): ...
x: TileTensor[DType.float32, ...]        # `...` always goes last
def g[dtype: DType, //, filter_len: Int](...)   # infer-only params before `//`
```

Check the parameter's default before dropping it — removing an explicit
`[2]` where the default is `1` changes behavior.

A default is not a constraint: with `b_type: DType = a_type`, a check
comparing the two never fires for callers that omit `b_type`.

## 4. Tuples and variadic `Coord`, not `IndexList`/`Index`

Tuples convert implicitly to coordinate types. A single-element tuple is
`(x,)`.

```mojo
# WRONG                                   # CORRECT
Index(m, n)                               (m, n)
f(Index(tokens))                          f((tokens,))
t.load[width=1](IndexList[2](i, k))       t.load((i, k))
Coord(IndexList[2](a, b))                 Coord(a, b)
Coord(Idx[1], Idx[5])                     coord[1, 5]
t.load((Idx[5], Idx[3]))                  t.load(coord[5, 3])
```

If a mixed-type tuple fails to unify (e.g. `Tuple[Int, UInt32]`), cast the
odd element; don't fall back to `IndexList`. `Coord` accepts both
`Int32`- and `Int64`-backed index lists, so don't add canonicalizing
conversions. Prefer the `coord` comptime helper for completely static
shapes and sizes.

## 5. `comptime` for compile-time-known work

```mojo
# WRONG                                   # CORRECT
for i in range(CHUNK):                    comptime for i in range(CHUNK):
nested comptime for x3                    comptime for a, b, c in product(range(X), range(Y), range(Z)):  # std.itertools
for i in range(rank - 1, -1, -1):         comptime for i in reversed(range(rank)):
var scale = 1.0 / sqrt(Float32(D))        comptime scale = Float32(1.0) / sqrt(Float32(D))
if ctx.api() == "metal":  (runtime)       comptime if ...:   # target/vendor checks are always comptime
```

Don't unroll a loop that doesn't benefit — reviewers ask "do you need to
unroll here?". If a comptime `Float32` math call fails to fold on GPU,
file an issue rather than silently changing the dtype.

Keep comptime dimensions static all the way into layouts. Dispatch
heuristics read `static_shape`; passing a comptime value as a runtime
`Int` erases it and silently routes to a slower path.

```mojo
# WRONG: N is comptime, but c.static_shape[1] == -1
var c = TileTensor(ptr, row_major((M, N)))
# CORRECT
var c = TileTensor(ptr, row_major((M, Idx[N])))
```

## 6. Don't `rebind` unless the types already match

`rebind` is a compile-time promise that two types are identical after
elaboration, not a cast. Try deleting every `rebind` you add or touch; keep
it only if the code stops building without it. To reinterpret a pointer, use
`ptr.unsafe_bitcast[T]()`; to convert a `SIMD` value, use
`x.cast[DType]()`.

```mojo
# WRONG: rebinds a type back to itself
p.init_pointee_copy(rebind[StrideType](StrideType()))
# CORRECT
p.init_pointee_copy(StrideType())

# WRONG: bridges different TileTensor layouts; the stored Coord shapes the
# ABI, so row_major[4, 4]() and row_major((Idx[4], 4)) differ in size
rebind[TileTensor[dtype, LayoutB, origin]](tile_a)
# CORRECT: build the tile with the layout the consumer expects

# WRONG: spell out params, then rebind so inference works
# CORRECT: pass tuples/values directly and let parameters infer
```

- After a migration, check whether old rebinds are still needed.
- Legitimate: `rebind[IndexList[rank]](x.canonicalize())` (valid only
  after `canonicalize()`), `SIMD` width/dtype reconciliation inside
  epilogue lambdas, and return-boundary nominal mismatches.
- Use an existing coercion helper (e.g. `_coerce_dynamic[T]`) rather
  than repeating `rebind[T](Scalar[...](...))`.
- For a constraint the compiler can't prove, use a `where` clause or
  `comptime assert`, not a rebind that hides the mismatch. A `where`
  clause pushes the check up to the caller (caught by the LSP, but generic
  callers must restate it); a `comptime assert` pushes it down into the
  body (still caught, only at elaboration). Use `where` for APIs called
  mostly with concrete types (`UInt8(1) in (1, 2, 3)` requiring
  `T in Self.Ts`), and `comptime assert` for APIs meant for generic
  contexts, where a `where` clause would force every caller to repeat the
  constraint.
- More than one rebind per function, or one you can't justify in a
  sentence, is a red flag.

The same applies to origins: don't call `as_any_origin()` or use
`ImmutAnyOrigin` without need, and use `TileTensor.as_imm()` rather
than overriding `mut=False`.

## 7. Mark read-only data immutable

Inputs that are never written take immutable origins so the signature
shows which arguments are outputs. Question every `MutAnyOrigin`.

```mojo
# WRONG                                        # CORRECT
src: Pointer[Float32, MutAnyOrigin]            src: ImmPointer[Float32, _]
```

Loads from immutable data are invariant by default; don't pass
`invariant=True` or other no-op arguments. Prefer the named type prefixes
to describe the mutability of the underlying data such as `ImmPointer`,
`MutPointer`, `ImmTileTensor`, and `MutTileTensor` over the `mut=True`
or `mut=False` prefix.

## 8. TileTensor over LayoutTensor and raw pointers

New and migrated kernels use TileTensor, coordinate indexing, and slice
syntax instead of flat pointer offsets. Direct pointer access was a
workaround for LayoutTensor and NDBuffer performance cliffs (they carried
dead bytes); a fully static TileTensor has the same ABI as a pointer, so
index the tensor directly.

```mojo
# WRONG                                   # CORRECT
out.ptr[row * N + col] = v                out[row, col] = v
src.slice(...)                            src[:h, :w]
```

Don't index raw pointers or compute offsets from strides in kernels.
Some accelerators don't expose raw pointers (their tensors are backed by
a memref), so `ptr[offset]` and hand-written `i * stride + j` math don't
port to them. Let TileTensor own the layout: index it by coordinate, and
use `tile`, `vectorize`, and `distribute` for sub-views instead of
pointer arithmetic.

```mojo
# WRONG
var p = t.unsafe_ptr()
p[row * row_stride + col] = v
var tile_ptr = p + (tm * BM) * N + tn * BN
# CORRECT
t[row, col] = v
var tile = t.tile[BM, BN](tm, tn)
```

Don't use `TileTensor.ptr`. When you do need the pointer, call
`TileTensor.unsafe_ptr()`, and only when you know the tensor is backed by
a pointer (not, for example, a memref).

Remaining `raw_load`/`raw_store` uses need a tracking issue.

## 9. Simplify control flow; remove duplication

- Return early instead of nesting `if/else`; no `else` after `return`.
- Iterate collections directly (`for d in dtypes:`), not by index.
- Start accumulators at their identity (`accum = 0`) instead of peeling
  the first iteration.
- Derive lists from existing ones (`float_dtypes.join(integer_dtypes)`).
- Lift repeated blocks into a helper or `comptime` alias.
- Return a `Tuple` from a helper rather than making callers do `+1, +2`
  offset arithmetic; pack scalar kernel args into one `InlineArray` and
  buffer instead of several buffers or fake `Pointer`s.
- Prefer t-strings over string concatenation. For example:

  ```mojo
  # WRONG
  "expected " + String(expected_dtype) + " but got " + String(actual_dtype)
  # CORRECT
  t"expected {expected_dtype} but got {actual_dtype}"
  ```

## 10. Naming

- No vendor or device class in generic APIs (`CompilationTarget(MI355X)`,
  not a `gpu_`-prefixed variant).
- Acronyms in capitals: `TMem`, not `Tmem`.
- Don't shadow builtins (`min_value`, not `min`).
- Descriptive over cryptic (`tA`, `tR` get flagged).
- Always declare locals with `var`.

## 11. API design

- `Optional[T]` instead of sentinel defaults.
- One generic entry point (`is[NVIDIA_GPU]()`) over a growing set of
  `is_xyz()` methods.
- Take `StringSlice`/`Span`/`Iterable` rather than owning types.
- Never accept an argument and silently ignore it on some targets.
- Keep advanced, low-level functions out of module exports.
- For `call_location()`, use a thin `@inline(.always)` wrapper that forwards
  the location to a `@inline(.never)` `_impl`.
- Use named results with `rebind[type_of(out)](x)` instead of restating
  a long return type (only when the rebind is otherwise required).

## 12. Invariants and hardware facts

- Query constants (`WARP_SIZE`, device limits); don't hard-code `32` or
  `128`.
- State shape and config assumptions with `comptime assert` /
  `debug_assert`; don't remove existing bounds checks.
- In tests, use `NaN`/`Inf` sentinels, not arbitrary magic numbers, and
  use prime and tile-straddling lengths (`13`, `63`/`64`/`65`) so the
  masking and tail code runs.
- An out-of-bounds write fix needs a guard band after the output, filled
  with a value the kernel can't produce and checked afterwards; comparing
  in-range elements can't see the overrun.
- Device graph capture replays the host's choices, so host code must not
  pick a kernel or size a grid from per-request values or from device
  data copied back to the host, and must handle `M == 0` warmup launches.
- Declare `MAX_THREADS_PER_BLOCK_METADATA` to match the real block size.

## 13. Write portable kernels; gate only what's vendor-specific

Specializing a kernel for one vendor is fine when it uses that vendor's
hardware (TMA, WGMMA, tcgen05 on NVIDIA; MFMA, `lgkmcnt` waits on AMD;
inline asm). A kernel built only from generic primitives (`thread_idx`,
`block_dim`, `barrier`, `warp.sum`, `WARP_SIZE`, shared memory, TileTensor)
isn't specific to the vendor it was written on, so don't gate it or its
test to that vendor by default. This applies in both directions: a kernel
written and tuned on AMD is no more AMD-only than one written on NVIDIA is
NVIDIA-only.

```bzl
# WRONG: generic kernel, gated to whichever vendor it was written on
gpu_constraints = ["//:nvidia_gpu"],
gpu_constraints = ["//:amd_gpu"],
# CORRECT: any GPU (add an Apple `incompatible` select only if it fails there)
gpu_constraints = ["//:has_gpu"],
```

- Before adding a vendor gate, run the kernel on the other vendor too (an
  AMD MI300X/MI355X or NVIDIA H100/B200 dev box, or let the other CI lane
  run it). If it also has a host path, try dropping the GPU constraint
  entirely so it runs on every accelerator.
- When reviewing, if a vendor-gated kernel uses no vendor intrinsics, ask
  the author to try enabling it for all GPUs, or all accelerators.
- Keep the vendor-specific part small: put the fast path behind
  `comptime if is_nvidia_gpu()` / `is_amd_gpu()` and keep a generic
  fallback, rather than gating the whole kernel.
- Don't gate generic optimizations on hardware. Vectorized loads and
  stores, `vectorize`, unrolling, and wider per-thread tiles usually help
  every vendor, so apply them unconditionally and size them with
  `simd_width_of` instead of putting them behind `is_nvidia_gpu()` or
  `is_amd_gpu()`. Gate one only if a benchmark shows it regresses a vendor.
- If it really fails on another vendor, gate it with a `TODO(<ticket>)`
  naming the failure (a missing intrinsic or a numerics gap), not just
  "NVIDIA only" or "AMD only".
- Don't hard-code one vendor's facts (warp size 32 vs 64, shared memory
  size) in generic code; see §12.

## 14. Hygiene

- Every temporary hack or API abuse gets a `TODO(<ticket>)` linking a
  filed issue.
- Keep comments short: a line or two on *why*, never a paragraph. Don't
  explain what the code plainly shows; delete comments that add nothing
  (`# (TileTensor-based)`).
- Don't commit agent status or design notes next to source.
- Changelog entries describe the end state, not intermediate steps.
