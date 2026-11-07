# ABOUTME: NetworkVisualizer: the net.visualize component drawing a network's
# ABOUTME: species and reactions as a layered graph or a stoichiometry matrix
"""The ``net.visualize`` component."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

import numpy as np

from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ...core.reaction._evaluate import check_missing
from ._common import finish, text_width
from ._graph import draw_graph, legend_widths
from ._layout import graph_layout, matrix_layout, species_layout
from ._species import draw_species, involved, species_legend_widths
from ._stoichiometry import draw_stoichiometry, stoichiometry_legend_widths
from ._themes import NetworkTheme, resolve_theme
from ._widths import edge_widths, width_key

if TYPE_CHECKING:
    from ...core.network.network import Network
    from ...core.reaction._evaluate import FunctionValue

#: Default ``(min, max)`` width of rate-scaled edges, points.
ARROW_WIDTH = (0.4, 3.5)


class NetworkVisualizer:
    """Figures of a network's species and reactions (``net.visualize``).

    Both views are static matplotlib figures with fixed spacing: the figure
    grows with the network, so every label stays readable. Save large networks
    as PDF or SVG.

    Parameters
    ----------
    network : Network
        The network to draw.

    Attributes
    ----------
    network : Network
        The network drawn.
    """

    def __init__(self, network: Network):
        """Hold *network*; nothing is drawn until a view is asked for.

        Parameters
        ----------
        network : Network
            The network to draw.
        """
        self.network = network

    def _rate_widths(
        self,
        reactions: Sequence[int],
        rate_values: Mapping[str, Any] | None,
        rate_functions: Mapping[str, FunctionValue] | None,
        arrow_width: tuple[float, float],
    ) -> tuple[dict[int, float] | None, list[tuple[float, float]]]:
        """Edge width per reaction and the legend key; ``(None, [])`` unscaled."""
        if rate_values is None:
            if rate_functions is not None:
                raise ValueError("rate_functions needs rate_values")

            return None, []

        lo, hi = arrow_width
        if not 0 < lo <= hi:
            raise ValueError(
                f"arrow_width must be (min, max) with 0 < min <= max; got {arrow_width}"
            )

        drawn = [self.network.reactions[i] for i in reactions]
        check_missing(
            [r.rate_coefficient for r in drawn], rate_values, rate_functions or {}
        )
        k = np.array(
            [r.evaluate_rate_coefficient(rate_values, rate_functions) for r in drawn]
        )
        widths = edge_widths(k, lo, hi)

        return dict(zip(reactions, widths.tolist())), width_key(k, lo, hi)

    def graph(
        self,
        *,
        theme: str | NetworkTheme = "vivid",
        save: str | Path | None = None,
        show: bool = True,
        close: bool = True,
        tight: bool = False,
        rate_values: Mapping[str, Any] | None = None,
        rate_functions: Mapping[str, FunctionValue] | None = None,
        arrow_width: tuple[float, float] = ARROW_WIDTH,
    ) -> tuple[Figure, Axes]:
        """Draw the species and reactions as a layered graph.

        Species sit in columns by atom count; each reaction is a square, coloured
        by its type and labelled with its position in ``net.reactions``, with its
        equation above it. Grey lines lead from reactants into the reaction and
        arrows lead out to its products.

        Parameters
        ----------
        theme : str or NetworkTheme, optional
            Colour theme: a name in :data:`~jaff.plotting.NETWORK_THEMES` or a
            theme object, by default ``"vivid"``.
        save : str or Path, optional
            File to save the figure to; the format follows the suffix.
        show : bool, optional
            Show the figure (skipped on non-interactive backends), by default
            ``True``.
        close : bool, optional
            Close the figure in pyplot afterwards, by default ``True``; the
            returned figure still works.
        tight : bool, optional
            Crop the figure to its drawn content, with no padding around it, by
            default ``False``.
        rate_values : mapping, optional
            Scale each reaction's lines by its rate coefficient, evaluated for
            these values (``{"tgas": 100.0, "av": 1.0, ...}``; densities as one
            array under ``"nden"``), see
            :meth:`~jaff.core.reaction.Reaction.evaluate_rate_coefficient`.
            Widths go logarithmically from the smallest rate drawn to the
            largest; rates more than 12 decades below the largest get the
            minimum width. The legend shows up to three reference widths. By
            default the lines have fixed widths.
        rate_functions : mapping, optional
            With *rate_values*: a callable or constant per unknown function
            the rate coefficients use (e.g. a shielding-table lookup).
        arrow_width : tuple of float, optional
            ``(min, max)`` line width in points of rate-scaled lines, by
            default ``(0.4, 3.5)``.

        Returns
        -------
        tuple of (Figure, Axes)
            The figure and its single axes (inches, y pointing down).

        Raises
        ------
        ValueError
            If *theme* names no registered theme, a rate value is missing, or
            *arrow_width* is not ``(min, max)`` with ``0 < min <= max``.
        """
        th = resolve_theme(theme)
        everything = range(len(self.network.reactions))
        widths, key = self._rate_widths(
            everything, rate_values, rate_functions, arrow_width
        )
        legend = legend_widths(self.network, th)
        layout = graph_layout(self.network, text_width, legend, key_row=bool(key))
        fig, ax = draw_graph(self.network, layout, th, widths, key)

        return finish(fig, ax, theme=th, save=save, show=show, close=close, tight=tight)

    def stoichiometry(
        self,
        *,
        theme: str | NetworkTheme = "vivid",
        save: str | Path | None = None,
        show: bool = True,
        close: bool = True,
        tight: bool = False,
    ) -> tuple[Figure, Axes]:
        """Draw the net stoichiometry matrix: products minus reactants.

        One row per reaction (its equation, then its position in
        ``net.reactions`` in a circle coloured by type) and one column per
        species, grouped by atom count. A cell shows how many of the species the
        reaction consumes or produces; the theme's colours tell which.

        Parameters
        ----------
        theme : str or NetworkTheme, optional
            Colour theme: a name in :data:`~jaff.plotting.NETWORK_THEMES` or a
            theme object, by default ``"vivid"``.
        save : str or Path, optional
            File to save the figure to; the format follows the suffix.
        show : bool, optional
            Show the figure (skipped on non-interactive backends), by default
            ``True``.
        close : bool, optional
            Close the figure in pyplot afterwards, by default ``True``; the
            returned figure still works.
        tight : bool, optional
            Crop the figure to its drawn content, with no padding around it, by
            default ``False``.

        Returns
        -------
        tuple of (Figure, Axes)
            The figure and its single axes (inches, y pointing down).

        Raises
        ------
        ValueError
            If *theme* names no registered theme.
        """
        th = resolve_theme(theme)
        widths = stoichiometry_legend_widths(self.network, th)
        layout = matrix_layout(self.network, text_width, widths)
        fig, ax = draw_stoichiometry(self.network, layout, th)

        return finish(fig, ax, theme=th, save=save, show=show, close=close, tight=tight)

    def species(
        self,
        name: str,
        *,
        theme: str | NetworkTheme = "vivid",
        save: str | Path | None = None,
        show: bool = True,
        close: bool = True,
        tight: bool = False,
        rate_values: Mapping[str, Any] | None = None,
        rate_functions: Mapping[str, FunctionValue] | None = None,
        arrow_width: tuple[float, float] = ARROW_WIDTH,
    ) -> tuple[Figure, Axes]:
        """Draw the reactions that form and destroy one species.

        The species sits in the middle. On its left are the reactions forming
        it, each with its other reactants further left; on its right the
        reactions destroying it, each with its other products further right.
        Reactions are numbered squares coloured by type with their equation
        above, as in :meth:`graph`. A reaction with the species on both sides
        and no net change is left out.

        Parameters
        ----------
        name : str
            Name of the species, as in ``net.species[name]``.
        theme : str or NetworkTheme, optional
            Colour theme: a name in :data:`~jaff.plotting.NETWORK_THEMES` or a
            theme object, by default ``"vivid"``.
        save : str or Path, optional
            File to save the figure to; the format follows the suffix.
        show : bool, optional
            Show the figure (skipped on non-interactive backends), by default
            ``True``.
        close : bool, optional
            Close the figure in pyplot afterwards, by default ``True``; the
            returned figure still works.
        tight : bool, optional
            Crop the figure to its drawn content, with no padding around it, by
            default ``False``.
        rate_values : mapping, optional
            Scale each reaction's lines by its rate coefficient, evaluated for
            these values (``{"tgas": 100.0, "av": 1.0, ...}``; densities as one
            array under ``"nden"``), see
            :meth:`~jaff.core.reaction.Reaction.evaluate_rate_coefficient`.
            Widths go logarithmically from the smallest rate drawn to the
            largest; rates more than 12 decades below the largest get the
            minimum width. The legend shows up to three reference widths. By
            default the lines have fixed widths.
        rate_functions : mapping, optional
            With *rate_values*: a callable or constant per unknown function
            the rate coefficients use (e.g. a shielding-table lookup).
        arrow_width : tuple of float, optional
            ``(min, max)`` line width in points of rate-scaled lines, by
            default ``(0.4, 3.5)``.

        Returns
        -------
        tuple of (Figure, Axes)
            The figure and its single axes (inches, y pointing down).

        Raises
        ------
        KeyError
            If the network has no species *name*.
        ValueError
            If *theme* names no registered theme, a rate value is missing, or
            *arrow_width* is not ``(min, max)`` with ``0 < min <= max``.
        """
        th = resolve_theme(theme)
        specie = self.network.species[name]
        widths, key = self._rate_widths(
            involved(self.network, specie), rate_values, rate_functions, arrow_width
        )
        legend = species_legend_widths(self.network, specie, th)
        layout = species_layout(
            self.network, specie, text_width, legend, key_row=bool(key)
        )
        fig, ax = draw_species(self.network, layout, th, widths, key)

        return finish(fig, ax, theme=th, save=save, show=show, close=close, tight=tight)
