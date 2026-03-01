"""
JAFF Code Generator CLI Interface.

This module provides the command-line interface for the JAFF (Just Another File Format)
code generator. It processes template files containing JAFF directives and generates
code for chemical reaction networks in various programming languages (C, C++, Fortran,
Python, Rust, Julia, R).

Usage:
    python -m jaff.generate --network <network_file> [--outdir <dir>] [--indir <dir>] [--files <file1> <file2> ...]

Examples:
    # Generate code from a specific template directory
    python -m jaff.generate --network networks/react_COthin --indir templates/ --outdir output/

    # Use a predefined template collection
    python -m jaff.generate --network networks/react_COthin --template my_template --outdir output/

    # Process specific files with a default language
    python -m jaff.generate --network networks/test.dat --files file1.txt file2.txt --lang rust
"""

import argparse
import warnings
from pathlib import Path

from jaff import Codegen as cg
from jaff import Network
from jaff.file_parser import Fileparser


def main() -> None:
    """
    Main entry point for the JAFF code generator CLI.

    Parses command-line arguments, validates input files and directories, and processes
    template files to generate code based on the specified chemical reaction network.

    Command-line Arguments:
        --network: Path to the chemical reaction network file (required)
        --outdir: Output directory for generated files (optional, defaults to current directory)
        --indir: Input directory containing template files to process (optional)
        --files: Individual template files to process (optional)
        --template: Name of a predefined template directory to use (optional)
        --lang: Default programming language for files without language detection (optional)
    """
    # Set up argument parser with comprehensive help text
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog="jaff.generate",
        description="Generate code for chemical reaction networks in multiple programming languages.",
        epilog="""
Examples:
  # Generate from a template directory
  python -m jaff.generate --network networks/react_COthin --indir templates/ --outdir output/

  # Use a predefined template collection
  python -m jaff.generate --network networks/react_COthin --template chemistry_solver --outdir output/

  # Process specific files with Rust
  python -m jaff.generate --network networks/test.dat --files rates.txt odes.txt --lang rust --outdir output/

  # Combine template and custom files
  python -m jaff.generate --network networks/test.dat --template base --files custom.cpp --outdir output/

Supported Languages:
  c, cxx (c++, cpp), fortran (f90), python (py), rust (rs), julia (jl), r

For more information, visit: https://github.com/tgrassi/jaff
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        "--network",
        required=True,
        metavar="FILE",
        help="Path to chemical reaction network file (required)",
    )

    # Output options
    parser.add_argument(
        "--outdir",
        metavar="DIR",
        help="Output directory for generated files (default: jaff/generated)",
    )

    # Input source options (mutually compatible)
    input_group = parser.add_argument_group("Input sources (can combine multiple)")
    input_group.add_argument(
        "--indir", metavar="DIR", help="Directory containing template files to process"
    )
    input_group.add_argument(
        "--files",
        nargs="+",
        metavar="FILE",
        help="Individual template file(s) to process",
    )
    input_group.add_argument(
        "--template",
        metavar="NAME",
        help="Name of predefined template collection in jaff/templates/generator/",
    )

    # Code generation options
    gen_group = parser.add_argument_group("Code generation options")
    gen_group.add_argument(
        "--lang",
        metavar="LANGUAGE",
        choices=["c", "cxx", "fortran", "python", "rust", "julia", "r"],
        help="Default programming language for unsupported files (choices: %(choices)s)",
    )
    args: argparse.Namespace = parser.parse_args()

    # Extract command-line arguments
    output_dir: str | None = args.outdir
    input_dir: str | None = args.indir
    input_files: list[str] | None = args.files
    network_file: str | None = args.network
    default_lang: str | None = args.lang
    template: str | None = args.template

    # List to collect all files to process
    files: list[Path] = []

    # Locate JAFF package directory and built-in template directory
    # Templates are stored in jaff/templates/generator/
    jaff_dir: Path = Path(__file__).parent
    generator_template_dir: Path = jaff_dir / "templates" / "generator"
    preprocessor_template_dir: Path = jaff_dir / "templates" / "preprocessor"

    # Validate network file is provided
    if network_file is None:
        raise RuntimeError("No network file supplied. Please enter a network file")

    # Resolve and validate network file path
    netfile: Path = Path(network_file).resolve()
    if not netfile.exists():
        raise FileNotFoundError(f"Unable to find network file: {netfile}")

    if not netfile.is_file():
        raise FileNotFoundError(f"{netfile} is not a valid file")

    # Handle output directory
    if output_dir is None:
        warnings.warn(
            "\n\nNo output directory has been supplied.\n"
            f"Files will be generated at {Path.cwd()}"
        )

    outdir: Path = (
        Path(output_dir).resolve() if output_dir is not None else jaff_dir / "generated"
    )
    if not outdir.exists():
        # Create output directory if it doesn't exist
        Path.mkdir(outdir)

    if not outdir.is_dir():
        raise NotADirectoryError(f"Output path is not a directory: {outdir}")

    # Handle predefined template directory if specified
    if template is not None:
        # Get list of available template directory names
        # Each subdirectory in templates/generator/ is a template collection
        generator_templates: list[str] = [
            file.name for file in generator_template_dir.iterdir() if file.is_dir()
        ]

        # Validate that the requested template exists
        if template not in generator_templates:
            raise ValueError(
                f"Invalid template name. Supported templates are {generator_templates}"
            )

        # Recursively collect all files from the template directory
        generator_template_path: Path = generator_template_dir / template
        preprocessor_template_path: Path = preprocessor_template_dir / template
        generator_files = [
            file for file in generator_template_path.rglob("*") if not file.is_dir()
        ]
        preprocesor_files = [
            file for file in preprocessor_template_path.rglob("*") if not file.is_dir()
        ]

        # Keep preproc files that don't have a corresponding generator file, otherwise use the generator file
        generator_file_names = [file.name for file in generator_files]
        for file in preprocesor_files:
            if file.name in generator_file_names:
                files.append(generator_files[generator_file_names.index(file.name)])
                continue
            files.append(file)

    # Collect files from input directory if specified
    if input_dir is not None:
        indir: Path = Path(input_dir).resolve()
        files.extend([f for f in indir.iterdir() if f.is_file()])

    # Collect individual files if specified
    if input_files is not None:
        for file in input_files:
            infile: Path = Path(file).resolve()

            if not infile.exists():
                raise FileNotFoundError(f"Invalid file path {file}")

            if not infile.is_file():
                raise FileNotFoundError(f"{file} is not a file")

            files.append(infile)

    # Ensure at least one input file was provided
    if not files:
        raise RuntimeError("No valid input file/folder/template have been supplied")

    # Ensure default language is supported by jaff code generation
    if default_lang and default_lang not in cg.get_language_tokens():
        raise ValueError(f"Unsupported language specified: {default_lang}")

    # Create a new network instance for each file
    net: Network = Network(str(netfile))

    # Process each template file
    for file in files:
        # Initialize file parser for this template
        fparser: Fileparser = Fileparser(net, file, default_lang)

        # Parse and generate code
        lines: str = fparser.parse_file()

        # Write generated code to output file
        outfile: Path = outdir / file.name
        with open(outfile, "w") as f:
            f.write(lines)

        print(f"{file.name} created at {outdir}")

    print(f"\nSuccessfully generated files\nGenerated files can be found at {outdir}")


if __name__ == "__main__":
    main()
