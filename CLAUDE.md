# dinara-align

Exact DNA alignment in Mojo: affine-gap Gotoh alignment on CPU and GPU, and a bit-parallel unit-cost
edit distance after A*PA2.

## Mojo skills

Modular's Mojo skills live in `.claude/skills` and correct pretrained Mojo, which is out of date.
A session-start hook (`.claude/settings.json`) reminds the agent to load them; follow it:

- `mojo-syntax` before writing, editing or reviewing any Mojo, with `references/idiomatic-mojo.md`.
- `mojo-gpu-fundamentals` before touching GPU code; `mojo-python-interop` before Mojo-Python interop.
- Load a skill again after a context compaction, which truncates loaded skills.

`pixi run sync-skills` refreshes them from upstream; `.claude/skills/UPSTREAM` records the commit.

## Working here

- Build and test: `pixi run test`; format: `pixi run format` (120 columns).
- Finished work merges straight into `main` and is pushed; no pull requests.
