---
tags:
    - License
icon: phosphor/copyright
---

# License

JAFF is released under the **MIT License**.

## MIT License

```text
Copyright 2025 Jaff Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## What This Means

The MIT License is a **permissive** open-source license that allows you to:

### Permissions

- **Commercial Use**: Use JAFF in commercial projects
- **Modification**: Modify the source code to suit your needs
- **Distribution**: Distribute JAFF and your modifications
- **Private Use**: Use JAFF for private/internal projects
- **Sublicense**: Include JAFF in projects with different licenses

### Conditions

- **License and Copyright Notice**: Include the MIT License and copyright notice in all copies or substantial portions of the software

### Limitations

- **Liability**: The authors are not liable for any damages or issues
- **Warranty**: The software is provided "as is" without warranty of any kind

## Using JAFF in Your Project

### Attribution

When using JAFF in your project, please include:

1. **The MIT License text** (shown above)
2. **Copyright notice**: "Copyright 2025 Jaff Contributors"

You can include this in:

- Your project's LICENSE file
- A THIRD_PARTY_LICENSES file
- Your software's about/credits section
- Your documentation

### Example Attribution

**In a LICENSE or THIRD_PARTY_LICENSES file:**

```text
This project uses JAFF (Just Another Fancy Format):

Copyright 2025 Jaff Contributors
Licensed under the MIT License
https://github.com/jaff-chemistry/jaff
```

**In source code:**

```cpp
/*
 * This code was generated using JAFF
 * Copyright 2025 Jaff Contributors
 * Licensed under the MIT License
 * https://github.com/jaff-chemistry/jaff
 */
```

## Generated Code

The chemistry-specific output JAFF produces from *your* network — the rate, ODE,
and Jacobian expressions — is your work, and you may license it however you choose,
including in proprietary software.

The caveat is **bundled template code**: JAFF ships solver templates (under
`templates/`) that are themselves MIT-licensed, and some generated files
incorporate that template code. Where your output contains JAFF template code, that
portion remains under the MIT License and its copyright/permission notice should be
retained.

**Practical guidance:**

- [x] The generated equations derived from your network are yours to license.
- [x] Retain the MIT notice for any bundled template code included in the output.
- [x] Attribution to JAFF for the equation output is appreciated but not required.

!!! note
    This is a plain-language summary, not legal advice. The [MIT License](#mit-license)
    text above governs; consult a professional for specific licensing questions.

## Contributing

By contributing to JAFF, you agree that your contributions will be licensed under the MIT License.

### Contributor Agreement

When you submit code to JAFF:

1. You retain copyright to your contributions
2. You grant the project a perpetual, worldwide, non-exclusive, royalty-free license
3. Your contributions will be distributed under the MIT License

See the [Contributing Guide](../development/contributing.md) for more details.

## Acknowledgments

JAFF is built on the shoulders of giants. We thank:

- The **NumPy**, **SciPy**, **SymPy**, and **Pandas** communities
- The **Python** core developers
- The **open-source** community at large

## See Also

- [Contributing Guide](../development/contributing.md) - How to contribute to JAFF
- [GitHub Repository](https://github.com/jaff-chemistry/jaff) - Source code and issues
