# Development

## Development Setup

### Clone and Install

```bash
git clone https://github.com/tgrassi/jaff.git
cd jaff
pip install -e ".[dev]"
```

### Development Dependencies

The development installation includes:

- **pytest**: Testing framework
- **pytest-cov**: Coverage reporting  
- **black**: Code formatting
- **ruff**: Fast Python linter

## Code Organization

### Source Structure

```
src/jaff/
├── __init__.py          # Package initialization
├── cli.py               # Command-line interface
├── network.py           # Main Network class
├── parsers.py           # Format-specific parsers
├── reaction.py          # Reaction class
├── species.py           # Species class
├── builder.py           # Code generation entry (templates/plugins)
├── preprocessor.py      # Template preprocessor
├── function_parser.py   # Auxiliary functions file parser
└── data/
    └── atom_mass.dat    # Atomic mass data
templates/
  ├── python_solve_ivp/
  ├── fortran_dlsodes/
  └── kokkos_ode/
plugins/
  ├── python_solve_ivp/
  ├── fortran_dlsodes/
  └── kokkos_ode/
```

### Key Components

- **Network**: Main class that coordinates loading, validation, and analysis
- **Species**: Represents individual chemical species with mass, charge, composition
- **Reaction**: Represents chemical reactions with rate expressions
- **Parsers**: Format-specific parsing functions for different network types
- **Builder/Plugins**: Code generation entry point and per-template plugins
- **Preprocessor**: Lightweight preprocessor to inject generated code into templates
- **Function Parser**: Parser for optional auxiliary function files

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=jaff

# Run specific test file
pytest tests/test_network.py
```

### Test Structure

Tests should be organized in the `tests/` directory:

```
tests/
├── __init__.py
├── test_network.py     # Network class tests
├── test_species.py     # Species class tests
├── test_reaction.py    # Reaction class tests
├── test_parsers.py     # Parser function tests
└── fixtures/           # Test data files
```

## Code Quality

### Formatting with Black

```bash
# Format all code
black src/jaff

# Check formatting without changes
black --check src/jaff
```

### Linting with Ruff

```bash
# Lint all code
ruff src/jaff

# Fix auto-fixable issues
ruff --fix src/jaff
```

## Adding New Format Support

To add support for a new network format:

### 1. Create Parser Function

Add a new parser function in `parsers.py`:

```python
def parse_myformat(line):
    """Parse a line from MyFormat network file.
    
    Args:
        line (str): Line from network file
        
    Returns:
        tuple: (reactants, products, tmin, tmax, rate)
    """
    # Parse the line and extract components
    reactants = [...]  # List of reactant names
    products = [...]   # List of product names  
    tmin = ...         # Minimum temp (or None)
    tmax = ...         # Maximum temp (or None)
    rate = ...         # Rate expression string
    
    return reactants, products, tmin, tmax, rate
```

### 2. Add Format Detection

Update the format detection logic in `network.py` (around line 170):

```python
# Example format detection logic (add to existing if/elif chain)
if "kida_marker" in srow:
    rr, pp, tmin, tmax, rate = parse_kida(srow)
elif "udfa_marker" in srow:
    rr, pp, tmin, tmax, rate = parse_udfa(srow)
elif "myformat_marker" in srow:
    rr, pp, tmin, tmax, rate = parse_myformat(srow)
else:
    # Default fallback
    continue
```

### 3. Add Tests

Create tests for the new parser:

```python
def test_parse_myformat():
    line = "example myformat line"
    reactants, products, tmin, tmax, rate = parse_myformat(line)
    
    assert reactants == ["H", "H"]
    assert products == ["H2"]
    # etc.
```

### 4. Update Documentation

- Add format description to `docs/formats.md`
- Include examples in `docs/examples.md`
- Update README.md

## Release Process

### Version Management

Update version in `pyproject.toml`:

```toml
[project]
version = "0.2.0"
```

### Testing Before Release

```bash
# Run full test suite
pytest --cov=jaff

# Test installation
pip install -e .
jaff --help

# Test documentation build
mkdocs build
```

### Building Documentation

```bash
# Install documentation dependencies
pip install mkdocs mkdocs-material mkdocstrings[python]

# Serve locally
mkdocs serve

# Build static site
mkdocs build
```

## Contributing Guidelines

### Pull Request Process

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make changes with tests
4. Run code quality checks: `black src/jaff && ruff src/jaff`
5. Run tests: `pytest`
6. Submit pull request

### Code Style

- Follow PEP 8 style guidelines
- Use Black for consistent formatting
- Add docstrings for all public functions
- Include type hints where appropriate
- Keep functions focused and small

### Commit Messages

Use clear, descriptive commit messages:

```
Add support for MyFormat network parser

- Implement parse_myformat() function
- Add format detection logic
- Include tests and documentation
- Update CLI help with new format info
```
