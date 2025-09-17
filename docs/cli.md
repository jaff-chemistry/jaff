# Command Line Interface

The `jaff` command provides a convenient way to work with chemical networks from the command line.

## Usage

```bash
jaff [OPTIONS] NETWORK_FILE
```

## Arguments

- `NETWORK_FILE`: Path to the chemical network file to load

## Options

### Basic Options

- `-h, --help`: Show help message and exit
- `-l LABEL, --label LABEL`: Set a custom label for the network (default: filename)
- `-e, --errors`: Exit immediately if validation errors are found

### Information Display

- `--list-species`: List all species in the network with their masses and charges
- `--list-reactions`: List all reactions in the network

### Validation (Not Yet Implemented)

- `--check-mass`: Check mass conservation in reactions
- `--check-charge`: Check charge conservation in reactions

!!! warning "CLI Validation Features"
    The `--check-mass` and `--check-charge` options are referenced in the CLI but not yet implemented in the Network class.

## Examples

### Basic Network Loading

Load a network and see the automatic validation output:

```bash
jaff networks/react_COthin
```

### List Network Contents

Display all species with their properties:

```bash
jaff networks/react_COthin --list-species
```

Display all reactions:

```bash
jaff networks/react_COthin --list-reactions  
```

### Error Handling

Exit immediately if validation errors are found:

```bash
jaff networks/react_COthin --errors
```

Set a custom network label:

```bash
jaff networks/react_COthin --label "My CO Network"
```

## Output

The CLI will display:

1. **Welcome message** with a randomly selected "fancy" word
2. **Network loading progress** with file information
3. **Validation results** including:
   - Sink/source species detection
   - Missing electron recombinations for ions
   - Isomer detection (species with same composition)
   - Duplicate reaction detection
4. **ODE generation status**
5. **Summary statistics** about variables found and reactions loaded

## Exit Codes

- `0`: Success
- `1`: Error occurred (file not found, parsing error, validation failure with `--errors`)