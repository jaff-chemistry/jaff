"""High-level build orchestrator for JAFF plugin-based code generation.

This module provides :class:`Builder`, a convenience wrapper that locates a
named plugin, resolves the template directory, and delegates the full
code-generation pipeline to the plugin's ``main()`` entry point.

Plugin discovery
----------------
Plugins live in ``jaff.plugins.<template_name>.plugin`` and must expose a
``main(network, *, path_template, path_build)`` function.  The *template_name*
corresponds to a sub-directory of ``TEMPLATES_DIR / "preprocessor"``
(``jaff/templates/preprocessor/``).
"""

import sys
from pathlib import Path

from jaff.config import TEMPLATES_DIR


class Builder:
    """Orchestrate plugin-based code generation for a chemical network.

    :class:`Builder` is the high-level entry point for users who want to
    generate a complete, ready-to-run solver project from a parsed network by
    selecting a named template.  Each template is a plugin module
    ``jaff.plugins.<template>.plugin`` paired with a directory of template
    files under ``jaff/templates/preprocessor/<template>/``; the plugin
    generates the code blocks with :class:`~jaff.codegen.Codegen` and writes
    the filled-in files with :class:`~jaff.codegen.Preprocessor`.

    The bundled templates are:

    - ``python_solve_ivp``: a Python solver driven by
      ``scipy.integrate.solve_ivp``.
    - ``fortran_dlsodes``: a Fortran solver using the ODEPACK DLSODES
      integrator, with a Makefile.
    - ``kokkos_ode``: a C++ ODE system for a header-only integrator library
      (VODE), with a CMake build.
    - ``microphysics``: network, rate and RHS sources for an AMReX
      Microphysics network.

    For finer-grained control over the generated code, write a ``$$JAFF``
    template and render it with :class:`~jaff.codegen.TemplateParser` (or the
    ``jaffgen`` command) instead.

    Parameters
    ----------
    network : Network
        Parsed chemical reaction network passed unchanged to the plugin.

    Attributes
    ----------
    network : Network
        The network instance provided at construction time.

    Examples
    --------
    >>> from pathlib import Path
    >>> from jaff import Network
    >>> from jaff.codegen import Builder
    >>> builder = Builder(Network("COthin"))
    >>> path = builder.build(template="python_solve_ivp", output_dir="output/")
    Building network with template: python_solve_ivp
    Network built successfully using template 'python_solve_ivp'.
    Output files are located in: output/
    >>> sorted(p.name for p in Path(path).iterdir())
    ['commons.py', 'fluxes.py', 'main.py', 'ode.py', 'rates.py']
    """

    def __init__(self, network) -> None:
        """Store the network for later use by :meth:`build`.

        Parameters
        ----------
        network : Network
            Parsed chemical reaction network.
        """
        self.network = network

    def build(
        self, template: str = "python_solve_ivp", output_dir: str | None = None
    ) -> str:
        """Generate a complete solver project using a named plugin template.

        Resolves the template directory, dynamically imports the plugin module
        ``jaff.plugins.<template>.plugin``, and calls its ``main()`` function
        with the network and resolved paths.  The plugin preprocesses its
        template files into the output directory (created if needed) and
        copies the template's remaining files there unchanged.  Progress
        messages are printed to stdout.

        If the plugin cannot be imported, all available template names are
        printed to stdout and the process exits with code ``1``.

        Parameters
        ----------
        template : str, optional
            Name of the plugin template to use: ``"python_solve_ivp"``,
            ``"fortran_dlsodes"``, ``"kokkos_ode"`` or ``"microphysics"``.
            Must match a sub-directory under ``jaff/templates/preprocessor/``
            and a ``jaff.plugins.<template>.plugin`` module.
            Default ``"python_solve_ivp"``.
        output_dir : str or None, optional
            Destination directory for the generated files.  When ``None`` the
            current working directory is used.

        Returns
        -------
        str
            Directory containing the generated files: *output_dir* as given,
            or the absolute path of the current working directory when
            *output_dir* is ``None``.

        Raises
        ------
        SystemExit
            If the named plugin module cannot be imported (exit code ``1``).
        """
        print("Building network with template:", template)

        # Resolve the template directory bundled with the jaff.codegen package
        templates_dir = TEMPLATES_DIR / "preprocessor"
        path_template = str(templates_dir / template)

        # Resolve the output directory (default: current working directory)
        if output_dir is None:
            path_build = str(Path.cwd())
        else:
            path_build = output_dir

        # Dynamically import the plugin module using its fully-qualified name
        try:
            module = __import__(f"jaff.plugins.{template}.plugin", fromlist=["main"])
        except ImportError:
            print(f"Error: Template '{template}' not found. Available templates are:")
            # List all sub-directories in the preprocessor templates folder
            for entry in sorted(templates_dir.iterdir()):
                print(entry.name)
            sys.exit(1)

        # Delegate all code generation and file writing to the plugin's main()
        module.main(self.network, path_template=path_template, path_build=path_build)

        print(f"Network built successfully using template '{template}'.")
        print(f"Output files are located in: {path_build}")

        return path_build
