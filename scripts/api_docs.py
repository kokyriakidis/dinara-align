#!/usr/bin/env python3
"""Writes docs/API.md, the reference of everything `dinara_align` exports, from its docstrings as
`mojo doc` reads them.

    pixi run api-docs

The names come from the package's `__init__.mojo`, in its order, each from the module that defines
it: functions with every overload's signature, types with their fields, constants and methods, and
constants alone. A docstring is the one source; nothing here is written by hand.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "API.md"


def exported() -> list:
    """`(module, name)` for every name `__init__.mojo` imports, in order."""
    text = (ROOT / "dinara_align" / "__init__.mojo").read_text()
    out = []
    for match in re.finditer(r"^from \.(\w+) import (\(([^)]*)\)|([^\n]+))$", text, re.MULTILINE):
        names = match.group(3) if match.group(3) is not None else match.group(4)
        for name in re.split(r"[,\s]+", names):
            if name:
                out.append((match.group(1), name))
    return out


def prose(entry: dict) -> str:
    """A declaration's docstring, its summary and description as paragraphs."""
    parts = [entry.get("summary", "").strip(), entry.get("description", "").strip()]
    return "\n\n".join(part for part in parts if part)


def function_block(entry: dict, level: str, owner: str = "") -> list:
    """A function's section at heading `level`: each overload's signature and docstring, a static
    method's name qualified by its type, `owner`."""
    lines = []
    for overload in entry["overloads"]:
        signature = overload["signature"]
        if owner:
            signature = signature.replace("def ", f"def {owner}.", 1) if overload.get("isStatic") else signature
        lines.append(f"```mojo\n{signature}\n```\n")
        text = prose(overload)
        if text:
            lines.append(text + "\n")
    return [f"{level} `{entry['name']}`\n"] + lines


def struct_block(entry: dict) -> list:
    """A type's section: its signature and docstring, a table of its fields, its constants, and its
    public methods with `__init__`."""
    lines = [f"### `{entry['name']}`\n", f"```mojo\n{entry['signature']}\n```\n"]
    text = prose(entry)
    if text:
        lines.append(text + "\n")
    if entry["fields"]:
        lines.append("| field | type | |\n| :-- | :-- | :-- |")
        for field in entry["fields"]:
            summary = field.get("summary", "").replace("\n", " ")
            lines.append(f"| `{field['name']}` | `{field['type']}` | {summary} |")
        lines.append("")
    for alias in entry.get("aliases", []):
        lines.append(f"- `{entry['name']}.{alias['name']}` = `{alias.get('value', '')}`: {alias.get('summary', '')}")
    if entry.get("aliases"):
        lines.append("")
    for method in entry["functions"]:
        if method["name"].startswith("_") and method["name"] != "__init__":
            continue
        lines += function_block(method, "####", entry["name"])
    return lines


def main() -> None:
    """Writes `OUT` from `mojo doc`'s JSON, the exported functions, types and constants in turn."""
    raw = subprocess.run(["mojo", "doc", str(ROOT / "dinara_align")], check=True, capture_output=True, text=True)
    package = json.loads(raw.stdout)["decl"]
    modules = {module["name"]: module for module in package["modules"]}
    lines = [
        "# API reference",
        "",
        "Everything `dinara_align` exports, from its docstrings; regenerate with `pixi run api-docs`.",
        "",
        prose(modules["__init__"]),
        "",
    ]
    sections = {"function": [], "struct": [], "alias": []}
    seen = set()
    for module_name, name in exported():
        module = modules.get(module_name)
        if module is None or name in seen:
            continue
        seen.add(name)
        for function in module["functions"]:
            if function["name"] == name:
                sections["function"] += function_block(function, "###")
        for struct in module["structs"]:
            if struct["name"] == name:
                sections["struct"] += struct_block(struct)
        for alias in module["aliases"]:
            if alias["name"] == name:
                sections["alias"].append(f"- `{name}` = `{alias.get('value', '')}`: {prose(alias)}")
    lines += ["## Functions", ""] + sections["function"]
    lines += ["## Types", ""] + sections["struct"]
    lines += ["## Constants", ""] + sections["alias"] + [""]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines))
    print(f"{OUT}: {len(seen)} names")


if __name__ == "__main__":
    sys.exit(main())
