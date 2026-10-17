---
tags:
    - Development
icon: phosphor/git-pull-request
---

# Contributing to JAFF

Thank you for your interest in contributing to JAFF! This guide covers the
workflow: forking, branching, opening a pull request, and getting it through CI.

All work lands on `main` through pull requests. You never push directly to
`main` — you push a branch to your fork and open a PR from it.

## Ways to Contribute

- **Report bugs** — open an issue describing the problem and how to reproduce it.
- **Suggest features** — open an issue to discuss before building anything large.
- **Fix issues / add features** — pick an issue, then send a PR.
- **Improve docs** — corrections and clarifications are always welcome.

For major changes, open an issue first so the direction can be agreed before you
invest time.

## 1. Fork and Clone

Fork the repository on GitHub, then clone **your fork**:

```bash
git clone https://github.com/YOUR_USERNAME/jaff.git
cd jaff
```

Add the upstream repository so you can keep your fork in sync:

```bash
git remote add upstream https://github.com/jaff-chemistry/jaff.git
git remote -v   # origin = your fork, upstream = jaff-chemistry/jaff
```

Set up the development environment (editable install with dev dependencies) by
following the [Installation guide](installation.md).

## 2. Create a Branch

Always branch off an up-to-date `main`. Never commit to `main` directly.

```bash
git checkout main
git pull upstream main          # sync with upstream
git checkout -b feature/short-description
```

### Branch Naming

Prefix the branch with its purpose so its intent is clear at a glance:

| Prefix      | Use for                           | Example                          |
| ----------- | --------------------------------- | -------------------------------- |
| `feature/`  | New functionality                 | `feature/gpu-codegen`            |
| `bug-fix/`  | Bug fixes                         | `bug-fix/fortran-index-offset`   |
| `docs/`     | Documentation changes             | `docs/rewrite-installation`      |
| `refactor/` | Restructuring, no behavior change | `refactor/split-template-engine` |
| `test/`     | Test additions or fixes           | `test/network-edge-cases`        |
| `chore/`    | Maintenance, tooling, deps        | `chore/bump-numpy`               |

## 3. Make Your Changes

Keep changes focused — one logical change per branch makes review faster.

### Commit Messages

Write clear messages explaining _what_ and _why_, not _how_:

```bash
# Good
git commit -m "feat: add support for GPU code generation"
git commit -m "fix: correct index offset in Fortran codegen"
git commit -m "docs: rewrite installation guide"

# Bad
git commit -m "fixed stuff"
git commit -m "WIP"
```

Suggested format:

```
type: short description (≤50 chars)

Longer explanation if needed. What changed and why,
wrapped at ~72 characters.

Fixes #123
```

Commit types: `feat`, `fix`, `docs`, `refactor`, `test`, `chore`.

### Tests and Code Style

Before opening a PR, make sure tests pass and the code is formatted and linted:

```bash
pytest          # run the test suite
ruff check .    # lint
ruff format .   # format
```

- See the [Testing Guide](testing.md) for running, writing, and covering tests.
- See the [Code Style Guide](code-style.md) for formatting, naming, type hints,
  and docstring conventions.

## 4. Open a Pull Request

Push your branch to your fork and open a PR against `jaff-chemistry/jaff:main`:

```bash
git push origin feature/short-description
```

GitHub will print a link to open the PR. In the PR description:

- Summarize what changed and why.
- Reference any related issues (`Fixes #123`).
- Note anything reviewers should pay attention to.

### Pre-submit Checklist

- [x] Branched off `main` with a correct prefix
- [x] `pytest` passes locally
- [x] `ruff check .` and `ruff format .` are clean
- [x] Documentation updated if behavior or interfaces changed
- [x] PR description explains the change and links related issues

## 5. Pass CI

Every pull request to `main` must pass all CI checks before it can be merged.
Two workflows run automatically:

| Workflow                 | What it does                                                                                   |
| ------------------------ | ---------------------------------------------------------------------------------------------- |
| **Tests**                | Runs `pytest` (with coverage) across Linux, macOS, and Windows on Python 3.11, 3.12, and 3.13. |
| **Deploy Documentation** | Builds the docs with `zensical build` to catch broken builds and links.                        |

The **Tests** run includes the [example notebooks](../getting-started/examples.md):
`tests/test_example_notebooks.py` executes each one headless and fails if it does
not exit 0. Tests that build generated code (the Fortran and C++ solvers, and the
notebooks that compile them) run only on Linux.

Because tests run on three operating systems and three Python versions, keep
code portable (use `pathlib` over hard-coded paths, avoid version-specific
syntax). If a check fails, open the failing job's logs from the PR's **Checks**
tab, fix locally, and push again — CI re-runs on every push.

## 6. Address Review Feedback

A maintainer will review your PR. To respond to comments, commit on the same
branch and push — the PR updates automatically and CI re-runs:

```bash
git add .
git commit -m "address review comments"
git push origin feature/short-description
```

Keep your branch current with `main` if it falls behind:

```bash
git fetch upstream
git rebase upstream/main
git push --force-with-lease origin feature/short-description
```

## Getting Help

- **GitHub Issues** — bug reports and feature requests.
- **GitHub Discussions** — questions about the codebase or usage.

## License

By contributing, you agree that your contributions are licensed under the
project's [MIT License](../about/license.md).

## See Also

- [Installation Guide](installation.md)
- [Testing Guide](testing.md)
- [Code Style Guide](code-style.md)
