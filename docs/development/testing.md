---
tags:
    - Development
    - Testing
icon: phosphor/bug
---

# Testing

## Overview

JAFF uses pytest for testing. This guide covers how to run, write, and organize tests.

## Running Tests

The sections below collect the pytest invocations you'll use day to day.

### Basic Commands

Run the whole suite, or narrow down to a file, a single test, or a name pattern.

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/test_network_parsers.py

# Run specific test
pytest tests/test_network_parsers.py::test_load_krome

# Run tests matching pattern
pytest -k "network"
```

### Coverage

Measure how much of `jaff` the suite exercises and generate a browsable report.

```bash
# Run with coverage
pytest --cov=jaff

# Generate HTML coverage report
pytest --cov=jaff --cov-report=html

# View report
open htmlcov/index.html
```

### Markers

Markers let you select subsets of tests by category.

<!-- prettier-ignore -->
!!! warning "Markers"
    Markers have not yet been added to the suite. The commands below show the intended usage once markers are defined (see [Test Markers](#test-markers)).

```bash
# Run only fast tests
pytest -m "not slow"

# Run only slow tests
pytest -m slow

# Run only unit tests
pytest -m unit
```

## Test Organization

How the `tests/` directory is laid out and how individual test files are structured.

### Directory Structure

Tests live in a flat `tests/` directory, one file per area, with sample network
files under `tests/fixtures/`:

```
tests/
├── __init__.py
├── test_network_initialization.py   # Network construction
├── test_network_parsers.py          # Multi-format parsing (KROME, KIDA, …)
├── test_network_validation.py       # Duplicate / sink / isomer checks
├── test_network_edge_cases.py       # Empty / malformed networks
├── test_network_json.py             # .jaff round-trip
├── test_sympy_json.py               # SymPy ↔ JSON encoding
├── test_cse_generation.py           # CSE in generated code
├── test_languages.py                # Per-language code generation
├── test_ode_and_jac.py              # ODE + Jacobian generation
├── test_repo_networks_json_roundtrip.py
└── fixtures/                        # Sample input files
    ├── sample_krome.dat
    ├── sample_kida.dat
    ├── sample_kida_valid.dat        # + sample_kida_valid.jfunc
    ├── sample_prizmo.dat
    ├── sample_uclchem.dat
    ├── sample_udfa.dat
    ├── empty_network.dat
    ├── malformed_network.dat
    ├── test_cse.dat
    ├── test_jac.dat
    └── test_jac_dedt.dat            # + test_jac_dedt.jfunc
```

### Test File Structure

Group related tests into `Test<Subject>` classes with descriptive method names.

```python
"""Tests for Network parsing."""

import pytest
from jaff import Network


class TestNetwork:
    """Network class tests."""

    def test_load_basic_network(self):
        """Test loading a basic network."""
        net = Network("tests/fixtures/sample_krome.dat")
        assert len(net.species) > 0
        assert len(net.reactions) > 0

    def test_load_krome_format(self):
        """Test loading KROME format."""
        net = Network("tests/fixtures/sample_krome.dat")
        assert net.label is not None


class TestNetworkValidation:
    """Network validation tests."""

    def test_duplicate_reactions(self):
        """Test detection of duplicate reactions."""
        # Test implementation
        pass
```

## Writing Tests

Patterns for the most common kinds of test, from a plain assertion to parametrized cases.

### Basic Test

The simplest test loads a fixture and asserts on the resulting object.

```python
def test_network_species_count():
    """Test species are loaded."""
    net = Network("tests/fixtures/sample_krome.dat")
    assert len(net.species) > 0
    assert net.species.count == len(net.species)
```

### Testing Exceptions

Use `pytest.raises` to assert that bad input fails the way you expect.

```python
from jaff.errors import ParserError


def test_missing_file_raises_error():
    """Test that a missing file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        Network("nonexistent.dat")

def test_malformed_network_raises_error():
    """Test that a malformed network raises a ParserError."""
    with pytest.raises(ParserError):
        Network("tests/fixtures/malformed_network.dat")
```

### Parametrized Tests

Run the same test body over many inputs with `@pytest.mark.parametrize`.

```python
@pytest.mark.parametrize("lang,expected_bracket", [
    ("c++", "["),
    ("c", "["),
    ("fortran", "("),
    ("python", "["),
])
def test_codegen_brackets(lang, expected_bracket):
    """Test bracket style for each language."""
    net = Network("tests/fixtures/sample_krome.dat")
    cg = Codegen(network=net, lang=lang)
    assert cg.lb == expected_bracket
```

### Parametrized with Multiple Arguments

Parametrize over several arguments at once to cover a matrix of cases.

```python
@pytest.mark.parametrize("filename,expected_species,expected_reactions", [
    ("sample_krome.dat", 10, 20),
    ("sample_kida.dat", 50, 100),
])
def test_network_sizes(filename, expected_species, expected_reactions):
    """Test networks load with expected sizes (counts are illustrative)."""
    net = Network(f"tests/fixtures/{filename}")
    assert len(net.species) == expected_species
    assert len(net.reactions) == expected_reactions
```

## Fixtures

Fixtures provide reusable setup. Define shared fixtures in a `conftest.py`, or
inline in the test module that uses them.

### Basic Fixtures

A fixture returns a prepared object that any test can request by name.

```python
# conftest.py
import pytest
from jaff import Network

@pytest.fixture
def small_network():
    """Small test network."""
    return Network("tests/fixtures/sample_krome.dat")

@pytest.fixture
def codegen_cpp(small_network):
    """C++ code generator for the small network."""
    from jaff import Codegen
    return Codegen(network=small_network, lang="c++")
```

### Using Fixtures

Request a fixture by adding its name as a test argument.

```python
def test_with_fixture(small_network):
    """Test using the network fixture."""
    assert len(small_network.species) > 0

def test_code_generation(codegen_cpp):
    """Test rate code generation."""
    rates = codegen_cpp.get_rates_str(use_cse=True)
    assert "k[0]" in rates
```

### Fixture Scope

Control how often a fixture is rebuilt with `scope` — reuse expensive setup
across a module, or get a fresh instance per test.

```python
@pytest.fixture(scope="module")
def expensive_network():
    """Module-scoped fixture for expensive setup (loaded once per module)."""
    return Network("tests/fixtures/sample_kida.dat")

@pytest.fixture(scope="function")
def temp_directory(tmp_path):
    """Function-scoped temporary directory (new for each test)."""
    return tmp_path
```

## Test Data

Where test inputs come from and how to create throwaway files.

### Creating Test Networks

Fixture network files are small plain-text reaction lists checked into `tests/fixtures/`.

```python
# tests/fixtures/sample_krome.dat
H + O -> OH, 1.2e-10 * (tgas/300)**0.5
H2 + O -> OH + H, 3.4e-11 * exp(-500/tgas)
C + O2 -> CO + O, 5.6e-12
```

### Using Temporary Files

Use pytest's `tmp_path` fixture for files a test writes and reads back.

```python
def test_network_save_load(tmp_path):
    """Test saving and loading a network."""
    # Create network
    net = Network("tests/fixtures/sample_krome.dat")

    # Save to a temporary .jaff file
    output_file = tmp_path / "test.jaff"
    net.to_jaff(output_file)

    # Load and verify
    net2 = Network(str(output_file))
    assert len(net2.species) == len(net.species)
```

## Mocking

Replace real dependencies with stand-ins to isolate the code under test.

### Mock External Dependencies

Patch I/O and other externals so tests don't touch the real filesystem.

```python
from unittest.mock import Mock, patch

def test_with_mock():
    """Test with a mocked open()."""
    with patch("builtins.open") as mock_open:
        mock_open.return_value.__enter__.return_value.readlines.return_value = [
            "H + O -> OH, 1.2e-10\n"
        ]
        # Test code
```

### Mock Network Loading

Build a `Mock(spec=Network)` when a test only needs a stand-in network.

```python
@pytest.fixture
def mock_network():
    """Mock network for testing."""
    net = Mock(spec=Network)
    net.species = [Mock(name="H"), Mock(name="O")]
    net.reactions = [Mock()]
    return net

def test_with_mock_network(mock_network):
    """Test using the mocked network."""
    assert len(mock_network.species) == 2
```

## Testing Best Practices

Habits that keep the suite readable and useful.

### 1. Test One Thing at a Time

Each test should verify a single behaviour so failures point to one cause.

```python
# Good
def test_network_loads_successfully():
    """Test network loads without error."""
    net = Network("tests/fixtures/sample_krome.dat")
    assert net is not None

def test_network_has_species():
    """Test network has species."""
    net = Network("tests/fixtures/sample_krome.dat")
    assert len(net.species) > 0

# Avoid
def test_everything():
    """Test everything at once."""
    net = Network("tests/fixtures/sample_krome.dat")
    assert net is not None
    assert len(net.species) > 0
    # ... many more assertions
```

### 2. Use Descriptive Names

Name the test after the behaviour it checks.

```python
# Good
def test_codegen_raises_error_for_unsupported_language():
    """Test that an unsupported language raises ValueError."""
    pass

# Avoid
def test_lang():
    pass
```

### 3. Test Edge Cases

Cover the boundaries — empty, single-element, and duplicate inputs.

```python
def test_empty_network():
    """Test handling of an empty network."""
    pass

def test_network_with_one_species():
    """Test a network with a single species."""
    pass

def test_network_with_duplicate_species():
    """Test handling of duplicate species."""
    pass
```

### 4. Test Error Conditions

Assert that invalid input raises the right exception.

```python
def test_file_not_found():
    """Test FileNotFoundError for a missing file."""
    with pytest.raises(FileNotFoundError):
        Network("nonexistent.dat")

def test_invalid_format():
    """Test ParserError for an invalid format."""
    from jaff.errors import ParserError
    with pytest.raises(ParserError):
        Network("tests/fixtures/malformed_network.dat")
```

## Integration Tests

Exercise several components together to confirm a full workflow holds up.

### Testing Complete Workflows

This test drives the whole path: load a network, generate code, write it out, verify.

```python
def test_complete_code_generation_workflow(tmp_path):
    """Test the complete workflow from network to code."""
    # Load network
    net = Network("tests/fixtures/sample_krome.dat")

    # Create code generator
    cg = Codegen(network=net, lang="c++")

    # Generate code
    rates = cg.get_rates_str(use_cse=True)
    odes = cg.get_ode_str(use_cse=True)

    # Save to file
    output_file = tmp_path / "chemistry.cpp"
    with open(output_file, "w") as f:
        f.write(rates)
        f.write("\n\n")
        f.write(odes)

    # Verify file exists and has content
    assert output_file.exists()
    content = output_file.read_text()
    assert "k[0]" in content
    assert "f[0]" in content
```

## Performance Tests

Guard against regressions in speed and memory. Mark these `slow` so they can be skipped.

### Timing Tests

Assert that an operation finishes within a time budget.

```python
import time

@pytest.mark.slow
def test_large_network_performance():
    """Test that a large network loads in reasonable time."""
    start = time.time()
    net = Network("tests/fixtures/sample_kida.dat")
    duration = time.time() - start

    assert duration < 30.0  # Should load in < 30 seconds
```

### Memory Tests

Track peak allocation with `tracemalloc` and assert an upper bound.

```python
import tracemalloc

@pytest.mark.slow
def test_memory_usage():
    """Test that memory usage stays reasonable."""
    tracemalloc.start()

    net = Network("tests/fixtures/sample_krome.dat")
    cg = Codegen(network=net, lang="c++")
    rates = cg.get_rates_str(use_cse=True)

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Peak memory should be reasonable
    assert peak < 100 * 1024 * 1024  # < 100 MB
```

## Continuous Integration

The **Tests** workflow runs `pytest` (with coverage) on every push and pull
request to `main`, across Linux, macOS, and Windows on Python 3.11, 3.12, and
3.13, and uploads coverage to Codecov. The suite also executes the
[example notebooks](../getting-started/examples.md) headless
(`tests/test_example_notebooks.py`) and builds and runs the solvers generated by
the Fortran and C++ templates (`tests/test_generated_solvers_run.py`); these build
tests run only on Linux. A second workflow builds the docs.

For the full list of CI workflows and what each one gates, see the
[Contributing Guide](contributing.md#5-pass-ci).

## Coverage Goals

Target coverage levels for the suite:

- **Overall coverage**: > 80%
- **Critical modules**: > 90%
- **New code**: 100% coverage required

### Check Coverage

Use `--cov-report=term-missing` to see exactly which lines are untested.

```bash
# Show untested lines in the terminal
pytest --cov=jaff --cov-report=term-missing

# Generate a browsable HTML report
pytest --cov=jaff --cov-report=html
open htmlcov/index.html
```

## Troubleshooting Tests

Flags and techniques for diagnosing failing tests.

### Test Failures

Increase verbosity, surface prints, or stop early to zero in on a failure.

```bash
# Run with detailed output
pytest -vv

# Show print statements
pytest -s

# Stop at first failure
pytest -x

# Run last failed tests
pytest --lf
```

### Debugging Tests

Drop into a debugger at the point of interest.

```python
def test_debug_example():
    """Test with a debugger breakpoint."""
    net = Network("tests/fixtures/sample_krome.dat")

    # Built-in breakpoint (Python 3.7+)
    breakpoint()

    assert len(net.species) > 0
```

## Test Markers

Markers tag tests so they can be selected or skipped as a group.

<!-- prettier-ignore -->
!!! warning "Markers"
    Markers are not yet defined in `pyproject.toml`. The snippets below show how they will be configured and used once added.

### Defining Markers

Register markers under `[tool.pytest.ini_options]` in `pyproject.toml`.

```toml
# pyproject.toml
[tool.pytest.ini_options]
markers = [
    "slow: marks tests as slow",
    "unit: unit tests",
    "integration: integration tests",
    "network: network-related tests",
]
```

### Using Markers

Apply markers with the `@pytest.mark.<name>` decorator; multiple markers stack.

```python
@pytest.mark.slow
def test_large_network():
    """Slow test."""
    pass

@pytest.mark.unit
def test_species_creation():
    """Unit test."""
    pass

@pytest.mark.integration
@pytest.mark.network
def test_full_workflow():
    """Integration test for the network workflow."""
    pass
```

## See Also

- [Contributing Guide](contributing.md)
- [Code Style Guide](code-style.md)
- [pytest Documentation](https://docs.pytest.org/)
