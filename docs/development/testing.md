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

Run the whole suite, or narrow down to a package, a file, a single test, or a name
pattern. Run pytest through `uv run --extra dev` so the dev tools are installed
(see the [Installation Guide](installation.md)).

```bash
# Run all tests
uv run --extra dev pytest

# Run with verbose output
uv run --extra dev pytest -v

# Run one package's tests
uv run --extra dev pytest tests/core/parsers

# Run specific test file
uv run --extra dev pytest tests/core/parsers/test_network_parsers.py

# Run specific test (Class::method)
uv run --extra dev pytest tests/core/parsers/test_network_parsers.py::TestNetworkParsers::test_krome_format_detection

# Run tests matching pattern
uv run --extra dev pytest -k "network"
```

### Coverage

Measure how much of `jaff` the suite exercises and generate a browsable report.

```bash
# Run with coverage
uv run --extra dev pytest --cov=jaff

# Generate HTML coverage report
uv run --extra dev pytest --cov=jaff --cov-report=html

# View report
open htmlcov/index.html
```

### Environment Variables

The suite defines no custom pytest markers. Optional behaviour is switched on with
environment variables instead:

| Variable                            | Effect                                                                                                   |
| ----------------------------------- | -------------------------------------------------------------------------------------------------------- |
| `JAFF_UPDATE_GOLDEN=1`              | Rewrite the committed generated code in `tests/golden/` from the current templates.                      |
| `JAFF_UPDATE_SNAPSHOTS=1`           | Record the snapshots (parser output in `tests/snapshots/parser_output/`, template-engine output in `tests/snapshots/template_engine/`) instead of comparing against them. |
| `JAFF_GOLDEN_SEED=<seed>`           | Fix the random input seed of the rendered-code tests, to reproduce a failure (the seed is in the failure message). |
| `JAFF_OFFLINE=1`                    | Do not contact the data server: no registry fetch and no downloads, so only already-installed data tables are available. |

After refreshing golden code or snapshots, review the diff before keeping it.

## Test Organization

How the `tests/` directory is laid out and how individual test files are structured.

### Directory Structure

Tests mirror the source tree: a test lives in the folder of the package whose
code it tests. Cross-cutting checks have their own folders, and test data sits
at the top level:

```text
tests/
├── conftest.py          # shared fixtures (make_network, quiet logging, …)
├── _paths.py            # TESTS_DIR, REPO_ROOT, FIXTURES_DIR, SNAPSHOTS_DIR, GOLDEN_DIR
├── helpers/             # shared helpers imported as tests.helpers.<module>: code-generation
│                        #   runners, build harnesses, test-only plugins, data builders
├── core/
│   ├── network/         # Network construction, validation, serialization, densities
│   ├── reaction/        # rate expressions, rate segments, temperature cutoffs
│   ├── species/         # Specie / Species
│   ├── elements/        # Elements
│   └── parsers/         # file formats (KIDA, KROME, PRIZMO, UCLCHEM, UDFA), .jfunc
├── physics/
│   ├── photo_reactions/ # radiation, cross sections, shielding, aux functions
│   ├── dust/
│   └── thermodynamics/  # EOS, internal energy, dE/dt and dT/dt
├── codegen/             # Codegen, languages, Jacobians
│   └── templating/      # the $$JAFF template engine
├── config/              # config loading, JSON schemas, reference page
├── cli/                 # jaffgen / jaffx
├── common/  types/  drivers/  utils/  plotting/
├── e2e/                 # golden code, generated solvers, bundled networks, snapshots
├── meta/                # public API, API docs, import side effects, notebooks
├── fixtures/            # sample input files
├── golden/              # committed generated code (JAFF_UPDATE_GOLDEN=1 to refresh)
└── snapshots/           # characterization snapshots
```

A module's config tests sit with the module (`physics/dust/test_dust_config.py`);
`config/` holds only the loading machinery. Read data through `tests._paths`
(`FIXTURES_DIR`, `SNAPSHOTS_DIR`, `GOLDEN_DIR`, `REPO_ROOT`, `TESTS_DIR`) or the
`fixtures_dir` fixture, never `Path(__file__)` or a path relative to the current
directory, so a test file can move — and pytest can run from anywhere — without
breaking.

### Test File Structure

Group related tests into `Test<Subject>` classes with descriptive method names.

```python
# ABOUTME: Tests for Network parsing
# ABOUTME: Loads the KROME sample fixture and checks the result

from jaff import Network
from tests._paths import FIXTURES_DIR


class TestNetwork:
    """Network class tests."""

    def test_load_basic_network(self):
        """Test loading a basic network."""
        net = Network(str(FIXTURES_DIR / "sample_krome.dat"))
        assert len(net.species) > 0
        assert len(net.reactions) > 0

    def test_label_is_file_stem(self):
        """The network label defaults to the file name without extension."""
        net = Network(str(FIXTURES_DIR / "sample_krome.dat"))
        assert net.label == "sample_krome"
```

## Writing Tests

Patterns for the most common kinds of test, from a plain assertion to parametrized cases.

### Basic Test

The simplest test requests a fixture and asserts on the resulting object.

```python
def test_network_species_count(sample_network):
    """Test species are loaded."""
    assert len(sample_network.species) > 0
    assert sample_network.species.count == len(sample_network.species)
```

### Testing Exceptions

Use `pytest.raises` to assert that bad input fails the way you expect.

```python
import pytest

from jaff import Network
from jaff.errors import ParserError
from tests._paths import FIXTURES_DIR


def test_missing_file_raises_error():
    """Test that a missing file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        Network("nonexistent.dat")


def test_malformed_network_raises_error():
    """Test that a malformed network raises a ParserError."""
    with pytest.raises(ParserError):
        Network(str(FIXTURES_DIR / "malformed_network.dat"))
```

### Parametrized Tests

Run the same test body over many inputs with `@pytest.mark.parametrize`.

```python
from jaff import Codegen


@pytest.mark.parametrize(
    "lang,expected_bracket",
    [
        ("c++", "["),
        ("c", "["),
        ("fortran", "("),
        ("python", "["),
    ],
)
def test_codegen_brackets(sample_network, lang, expected_bracket):
    """Test bracket style for each language."""
    cg = Codegen(network=sample_network, lang=lang)
    assert cg.lang.lb == expected_bracket
```

### Parametrized with Multiple Arguments

Parametrize over several arguments at once to cover a matrix of cases.

```python
@pytest.mark.parametrize(
    "filename,expected_species,expected_reactions",
    [
        ("sample_krome.dat", 9, 6),
        ("sample_kida.dat", 6, 7),
    ],
)
def test_network_sizes(filename, expected_species, expected_reactions):
    """Test networks load with the expected sizes."""
    net = Network(str(FIXTURES_DIR / filename))
    assert len(net.species) == expected_species
    assert len(net.reactions) == expected_reactions
```

## Fixtures

Fixtures provide reusable setup. Define shared fixtures in a `conftest.py`, or
inline in the test module that uses them.

### Shared Fixtures

`tests/conftest.py` provides these to every test:

| Fixture          | Scope    | What it gives you                                                                                  |
| ---------------- | -------- | -------------------------------------------------------------------------------------------------- |
| `_quiet`         | autouse  | Silences the JAFF banner (`print`) and logging, so test output stays clean.                         |
| `fixtures_dir`   | function | `tests/fixtures/` as a `Path` (the same as `tests._paths.FIXTURES_DIR`).                            |
| `make_network`   | function | Factory `make_network(lines, name=..., **kwargs)`: writes reaction lines to a file in `tmp_path` and loads it as a `Network`; `kwargs` go to `Network`. |
| `sample_network` | function | A `Network` loaded from `tests/fixtures/sample_kida.dat`.                                           |
| `snapshot_store` | function | Compares values against JSON snapshots in `tests/snapshots/parser_output/` (`JAFF_UPDATE_SNAPSHOTS=1` records them). |
| `rendered`       | session  | Renders the templates in `tests/fixtures/templates/` once per network and returns the output folder. |
| `render_seed`    | session  | Random input seed for rendered-code tests (`JAFF_GOLDEN_SEED` fixes it).                           |

`tests/plotting/conftest.py` adds, for the plotting tests only, an autouse fixture
that closes every Matplotlib figure after each test (and selects the `Agg` backend)
and `jaff_warnings`, a list that records the JAFF logger's warnings.

### Defining Fixtures

Define shared fixtures in a `conftest.py`, or inline in the test module that uses
them. A fixture can build on other fixtures:

```python
import pytest

from jaff import Codegen


@pytest.fixture
def codegen_cpp(sample_network):
    """C++ code generator for the sample network."""
    return Codegen(network=sample_network, lang="c++")
```

### Using Fixtures

Request a fixture by adding its name as a test argument.

```python
def test_with_fixture(sample_network):
    """Test using the network fixture."""
    assert len(sample_network.species) > 0


def test_inline_network(make_network):
    """Build a network from inline reaction lines."""
    net = make_network(["H + H -> H2 [10,1000] 1e-10"], name="my_net.dat")
    assert net.label == "my_net"


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
    return Network(str(FIXTURES_DIR / "sample_kida.dat"))


@pytest.fixture(scope="function")
def temp_directory(tmp_path):
    """Function-scoped temporary directory (new for each test)."""
    return tmp_path
```

## Test Data

Where test inputs come from and how to create throwaway files.

### Test Networks

Small network files in each supported format are checked into `tests/fixtures/`
(`sample_kida.dat`, `sample_krome.dat`, `sample_prizmo.dat`, `sample_uclchem.dat`,
`sample_udfa.dat`, plus edge cases such as `empty_network.dat` and
`malformed_network.dat`). For example, `sample_krome.dat` starts with:

```text
# Sample KROME network file
@format:idx,R,R,R,P,P,P,P,tmin,tmax,rate
@var:inv_tgas = 1.0/tgas
@var:t32 = tgas/3e2

# idx,R1,R2,R3,P1,P2,P3,P4,tmin,tmax,rate
1,H,H,H,H2,H,,,10,1000,5.0e-30
2,H2,PHOTON,,H,H,,,10,3000,1.0e-11
```

For a one-off network, prefer the `make_network` fixture over a new fixture file.

### Using Temporary Files

Use pytest's `tmp_path` fixture for files a test writes and reads back.

```python
def test_network_save(tmp_path, sample_network):
    """Test saving a network."""
    output_file = tmp_path / "test.jaff"
    sample_network.to_jaff(output_file)

    assert output_file.read_bytes()[:2] == b"\x1f\x8b"  # gzip
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
def test_network_loads_successfully(sample_network):
    """Test network loads without error."""
    assert sample_network is not None


def test_network_has_species(sample_network):
    """Test network has species."""
    assert len(sample_network.species) > 0


# Avoid
def test_everything():
    """Test everything at once."""
    net = Network(str(FIXTURES_DIR / "sample_krome.dat"))
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
        Network(str(FIXTURES_DIR / "malformed_network.dat"))
```

## Integration Tests

Exercise several components together to confirm a full workflow holds up.

### Testing Complete Workflows

This test drives the whole path: load a network, generate code, write it out, verify.

```python
def test_complete_code_generation_workflow(tmp_path):
    """Test the complete workflow from network to code."""
    # Load network
    net = Network(str(FIXTURES_DIR / "sample_kida.dat"))

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

## Continuous Integration

The **Tests** workflow installs with `uv sync --locked --extra dev` and runs `pytest` (with coverage) on every push and pull
request to `main`, across Linux, macOS, and Windows on Python 3.11, 3.12, and
3.13, and uploads coverage to Codecov. The suite also executes the
[example notebooks](../getting-started/examples.md) headless
(`tests/meta/test_example_notebooks.py`) and builds and runs the solvers generated by
the Fortran and C++ templates (`tests/e2e/test_generated_solvers_run.py`); these build
tests run only on Linux. A second workflow builds the docs.

For the full list of CI workflows and what each one gates, see the
[Contributing Guide](contributing.md#5-pass-ci).

## Coverage Goals

Target coverage levels for the suite:

- **Overall coverage**: > 80%
- **Critical modules**: > 90%
- **New code**: covered by tests written before the implementation (TDD)

### Check Coverage

Use `--cov-report=term-missing` to see exactly which lines are untested.

```bash
# Show untested lines in the terminal
uv run --extra dev pytest --cov=jaff --cov-report=term-missing

# Generate a browsable HTML report
uv run --extra dev pytest --cov=jaff --cov-report=html
open htmlcov/index.html
```

## Troubleshooting Tests

Flags and techniques for diagnosing failing tests.

### Test Failures

Increase verbosity, surface prints, or stop early to zero in on a failure.

```bash
# Run with detailed output
uv run --extra dev pytest -vv

# Show output (note: the autouse `_quiet` fixture replaces `print` with a no-op)
uv run --extra dev pytest -s

# Stop at first failure
uv run --extra dev pytest -x

# Run last failed tests
uv run --extra dev pytest --lf
```

### Debugging Tests

Drop into a debugger at the point of interest.

```python
def test_debug_example(sample_network):
    """Test with a debugger breakpoint."""
    net = sample_network

    # Built-in breakpoint (Python 3.7+)
    breakpoint()

    assert len(net.species) > 0
```

## See Also

- [Contributing Guide](contributing.md)
- [Code Style Guide](code-style.md)
- [pytest Documentation](https://docs.pytest.org/)
