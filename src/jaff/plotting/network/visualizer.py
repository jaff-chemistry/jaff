# ABOUTME: NetworkVisualizer: the net.visualize component drawing a network's
# ABOUTME: species and reactions as a layered graph or a stoichiometry matrix
"""The ``net.visualize`` component."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Sequence

from matplotlib.axes import Axes
from matplotlib.colors import Colormap, to_hex
from matplotlib.figure import Figure

from ...core.analysis._evaluate import evaluate
from ._common import finish, text_width
from ._graph import draw_graph, legend_widths
from ._layout import graph_layout, matrix_layout, species_layout
from ._species import draw_species, involved, species_legend_widths
from ._stoichiometry import draw_stoichiometry, stoichiometry_legend_widths
from ._themes import NetworkTheme, resolve_theme
from ._rates import RateColours, rate_colormap, rate_positions
from ._values import value_groups

if TYPE_CHECKING:
    from ...core.network.network import Network

#: What ``colour_by`` may name.
COLOUR_BY = ("rate_coefficient", "rate")


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

    def _rate_colours(
        self,
        reactions: Sequence[int],
        theme: NetworkTheme,
        colour_by: str | None,
        rate_cmap: str | Colormap | None,
    ) -> RateColours | None:
        """Rate colour of each drawn reaction, or ``None`` without *colour_by*."""
        if colour_by is None:
            if rate_cmap is not None:
                raise ValueError("rate_cmap needs colour_by")
            return None
        if colour_by not in COLOUR_BY:
            raise ValueError(
                f"colour_by must be one of {', '.join(COLOUR_BY)}, got {colour_by!r}"
            )
        if not self.network.analysis.enabled:
            raise ValueError(
                f"colour_by={colour_by!r} reads the values from analysis mode: "
                "call net.analysis.enable(...) first"
            )

        drawn = [self.network.reactions[i] for i in reactions]
        exprs = [r.rate if colour_by == "rate" else r.rate_coefficient for r in drawn]
        values = evaluate(self.network, exprs)
        positions, span = rate_positions(values)
        cmap = rate_colormap(theme, rate_cmap)
        colours = {i: to_hex(cmap(p)) for i, p in zip(reactions, positions)}

        return RateColours(
            colours,
            cmap,
            span if span is not None else (0.0, 0.0),
            colour_by,
            value_groups(self.network.analysis),
        )

    def graph(
        self,
        *,
        theme: str | NetworkTheme = "vivid",
        save: str | Path | None = None,
        show: bool = True,
        close: bool = True,
        tight: bool = False,
        colour_by: str | None = None,
        rate_cmap: str | Colormap | None = None,
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
        colour_by : {"rate_coefficient", "rate"}, optional
            Colour each reaction's marker and product arrows by its rate
            coefficient *k* or its rate 𝓡 (*k* times the reactant densities),
            evaluated with the values of ``net.analysis``, which must be
            enabled (see :class:`~jaff.core.analysis.Analysis`). Colours go
            logarithmically from the smallest value drawn to the largest along
            the theme's colormap; values more than 12 decades below the largest
            get the lowest colour. The reaction type is then shown by the marker
            shape alone, and a colour bar sits below the legend. By default
            reactions are coloured by type.
        rate_cmap : str or Colormap, optional
            With *colour_by*: a matplotlib colormap (or its name) to use, as
            given, instead of the theme's ``rates`` colormap.

        Returns
        -------
        tuple of (Figure, Axes)
            The figure and its single axes (inches, y pointing down).

        Raises
        ------
        ValueError
            If *theme* names no registered theme, *colour_by* is unknown or
            given with analysis mode off, a value is missing, or *rate_cmap* is
            given without *colour_by*.
        """
        th = resolve_theme(theme)
        everything = range(len(self.network.reactions))
        rated = self._rate_colours(everything, th, colour_by, rate_cmap)
        legend = legend_widths(self.network, th)
        layout = graph_layout(
            self.network,
            text_width,
            legend,
            colour_bar=rated is not None,
            values=rated.values if rated is not None else (),
        )
        fig, ax = draw_graph(self.network, layout, th, rated)

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
        colour_by: str | None = None,
        rate_cmap: str | Colormap | None = None,
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
        colour_by : {"rate_coefficient", "rate"}, optional
            Colour each reaction's marker and product arrows by its rate
            coefficient *k* or its rate 𝓡 (*k* times the reactant densities),
            evaluated with the values of ``net.analysis``, which must be
            enabled (see :class:`~jaff.core.analysis.Analysis`). Colours go
            logarithmically from the smallest value drawn to the largest along
            the theme's colormap; values more than 12 decades below the largest
            get the lowest colour. The reaction type is then shown by the marker
            shape alone, and a colour bar sits below the legend. By default
            reactions are coloured by type.
        rate_cmap : str or Colormap, optional
            With *colour_by*: a matplotlib colormap (or its name) to use, as
            given, instead of the theme's ``rates`` colormap.

        Returns
        -------
        tuple of (Figure, Axes)
            The figure and its single axes (inches, y pointing down).

        Raises
        ------
        KeyError
            If the network has no species *name*.
        ValueError
            If *theme* names no registered theme, *colour_by* is unknown or
            given with analysis mode off, a value is missing, or *rate_cmap* is
            given without *colour_by*.
        """
        th = resolve_theme(theme)
        specie = self.network.species[name]
        rated = self._rate_colours(
            involved(self.network, specie), th, colour_by, rate_cmap
        )
        legend = species_legend_widths(self.network, specie, th)
        layout = species_layout(
            self.network,
            specie,
            text_width,
            legend,
            colour_bar=rated is not None,
            values=rated.values if rated is not None else (),
        )
        fig, ax = draw_species(self.network, layout, th, rated)

        return finish(fig, ax, theme=th, save=save, show=show, close=close, tight=tight)
