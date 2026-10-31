# ABOUTME: NetworkVisualizer: the net.visualize component drawing a network's
# ABOUTME: species and reactions as a layered graph or a stoichiometry matrix
"""The ``net.visualize`` component."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ._common import finish, text_width
from ._graph import draw_graph, legend_widths
from ._layout import graph_layout, matrix_layout, species_layout
from ._species import draw_species, species_legend_widths
from ._stoichiometry import draw_stoichiometry, stoichiometry_legend_widths
from ._themes import NetworkTheme, resolve_theme

if TYPE_CHECKING:
    from ...core.network.network import Network


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

    def graph(
        self,
        *,
        theme: str | NetworkTheme = "vivid",
        save: str | Path | None = None,
        show: bool = True,
        close: bool = True,
        tight: bool = False,
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
        layout = graph_layout(self.network, text_width, legend_widths(self.network, th))
        fig, ax = draw_graph(self.network, layout, th)

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

        Returns
        -------
        tuple of (Figure, Axes)
            The figure and its single axes (inches, y pointing down).

        Raises
        ------
        KeyError
            If the network has no species *name*.
        ValueError
            If *theme* names no registered theme.
        """
        th = resolve_theme(theme)
        specie = self.network.species[name]
        widths = species_legend_widths(self.network, specie, th)
        layout = species_layout(self.network, specie, text_width, widths)
        fig, ax = draw_species(self.network, layout, th)

        return finish(fig, ax, theme=th, save=save, show=show, close=close, tight=tight)
