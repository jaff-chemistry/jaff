# Supported Network Formats

JAFF automatically detects and parses multiple astrochemical network formats. The detection is based on the file structure and specific format markers.

## Format Detection Logic

JAFF examines each line of the input file and applies format detection in this order:

1. **PRIZMO**: Lines containing `->` separator and `VARIABLES{}` blocks
2. **UDFA**: Lines containing `:` separator with specific column format
3. **KROME**: Comma-separated values with `@format:` header or multiple commas
4. **UCLCHEM**: Comma-separated values containing `,NAN,` marker
5. **KIDA**: Default format with fixed-width columns

## PRIZMO Format

Uses `->` as the reaction separator and supports custom variables.

### Features
- Reactions: `reactants -> products`
- Variable definitions in `VARIABLES{}` blocks
- Temperature limits: `TMIN`, `TMAX`

### Example
```
VARIABLES{
k1 = 1.0e-9
temp_limit = 300.0
}

H + H -> H2     k1    10.0    1000.0
```

## UDFA Format

Uses `:` as the separator with specific column structure.

### Features
- Reactions use `:` separator
- Fixed column format for rate parameters
- Temperature-dependent rate coefficients

### Example
```
H + H : H2 : 1.0e-9 : 0.0 : 0.0 : 10.0 : 1000.0
```

## KROME Format

Comma-separated values with optional format specification.

### Features
- Header: `@format:idx,R,R,R,P,P,P,P,tmin,tmax,rate`
- Variables: `@var:variable_name=expression`
- Rate expressions support mathematical functions
- Multiple reactants and products per reaction

### Example
```
@format:idx,R,R,R,P,P,P,P,tmin,tmax,rate
@var:k_formation=1.0e-9*(tgas/300.0)**0.5

1,H,H,,H2,,,10.0,1000.0,k_formation
```

## UCLCHEM Format

Comma-separated format with `,NAN,` marker for identification.

### Features
- Contains `,NAN,` somewhere in reaction lines
- Similar structure to KROME format
- Under construction - limited support

### Example
```
H,H,,H2,,,NAN,1.0e-9,10.0,1000.0
```

## KIDA Format (Default)

Fixed-width column format used as fallback.

### Features
- Fixed column positions for reaction components
- Whitespace-separated fields
- Rate coefficients with α, β, γ parameters: `k = α(T/300)^β exp(-γ/T)`

### Example
```
   1 H        H                   H2                    1.00E-09  0.50  0.00    10.0  1000.0
```

## Rate Expressions

All formats support rate expressions that can contain:

### Mathematical Functions
- Basic arithmetic: `+`, `-`, `*`, `/`, `**` (power)
- Functions: `exp()`, `sqrt()`, `log()`, `abs()`
- Constants: numerical values, `pi`, `e`

### Primitive Variables
- `tgas`: Gas temperature (K)
- `av`: Visual extinction (Draine units)  
- `crate`: Cosmic ray ionization rate (s⁻¹)
- `ntot`: Total number density (cm⁻³)
- `hnuclei`: H nuclei number density (cm⁻³)
- `d2g`: Dust-to-gas mass ratio

### Temperature Limits
- `tmin`: Minimum temperature for rate validity
- `tmax`: Maximum temperature for rate validity
- Applied as: `rate = rate.subs(tgas, max(tgas, tmin))` and `rate.subs(tgas, min(tgas, tmax))`

### Custom Variables
Each format allows definition of custom variables that are substituted into rate expressions before evaluation.

## Photo-chemistry

Special handling for photochemical reactions:

```
rate_expression = photo(reaction_index, cross_section, branching_ratio)
```

Photo-reactions are converted to function calls that can be evaluated with external photochemical rate data.