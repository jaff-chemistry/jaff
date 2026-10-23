---
tags:
    - User-guide
    - Code-generation
icon: phosphor/code-block
---

# Code Generation

JAFF can turn any loaded network into ready-to-compile/run source files in C, C++, Fortran, Python, Rust, Julia, or R. The generation pipeline is template-driven: you write ordinary source files that contain special `$$JAFF` directives, and JAFF expands them into network-specific code.

<div class="grid cards" markdown>

- :phosphor-textbox:{ .sm .middle } **Template Syntax**

    ***

    Learn the `$$JAFF` directive language: `SUB`, `REPEAT`, `REDUCE`, `GET`, `HAS`, `END`, option blocks with `REPLACE`, reusable `DEFINE`/`USE` sets, and the `LANG` and `COMMENT` file settings.

    [:octicons-arrow-right-24: Template Syntax](template-syntax.md)

- :phosphor-gear-six:{ .sm .middle } **Configuration File**

    ***

    Configure `jaffgen` runs declaratively with a `jaffgen.toml` file: set network paths, template names, output directories, radiation bands, and data-table conversions.

    [:octicons-arrow-right-24: Configuration File](../configuration/index.md)

- :phosphor-terminal-window:{ .sm .middle } **jaffgen CLI**

    ***

    Run the code-generation pipeline from the command line. Combine predefined templates, custom input directories, and individual files in a single invocation.

    [:octicons-arrow-right-24: jaffgen CLI](jaffgen.md)

- :phosphor-chart-line:{ .sm .middle } **Table Interpolation**

    ***

    Use precomputed lookup tables in rates and cooling: how `interp` functions in `.jfunc` are preserved through codegen, how their Jacobian derivatives become `_partial_N` calls, and how `[[table]]` emits the data.

    [:octicons-arrow-right-24: Table Interpolation](table-interpolation.md)

</div>
