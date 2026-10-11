"""Publication-quality plotting for JAFF, built on the seaborn objects API.

:class:`Plotter` renders tidy long DataFrames of labelled curves onto
caller-supplied or freshly created matplotlib axes via ``Plot.on(ax)``.  The
generic entry point is :meth:`Plotter.render_series`; the higher-level free
functions :func:`jaff.plotting.plot_rates` and :func:`jaff.plotting.plot_xsecs`
build frames and delegate to it.  Band bars and shaded fills are drawn with
matplotlib directly (the objects ``Area`` mark does not render reliably onto an
existing axes).

By default the house theme is applied *scoped* -- only while a figure is being
drawn -- so instantiating :class:`Plotter` never mutates global matplotlib
state.  Pass ``global_theme=True`` (or call
:func:`jaff.plotting.apply_global_theme`) to opt into a sticky session theme.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
import seaborn.objects as so

from . import _frames, _units, _xsec
from ._frames import _sampled_xsecs
from ._theme import (
    LOGO_PALETTE,
    apply_global_theme,
    despine,
    theme_context,
    theme_rc,
)

if TYPE_CHECKING:
    import pandas as pd

    from ..physics.photo_reactions._photochemistry import XsecsProps


class Plotter:
    """Publication-quality plotter using the seaborn objects interface.

    Parameters
    ----------
    palette : list[str] or None, optional
        Colour cycle for curves.  Defaults to the JAFF brand palette
        (:data:`jaff.plotting._theme.LOGO_PALETTE`).  Pass
        :data:`jaff.plotting._theme.MUTED_PALETTE` or
        :data:`jaff.plotting._theme.DEEP_PALETTE` for the seaborn cycles.
    global_theme : bool, optional
        If ``True``, apply the house theme globally and persistently on
        construction (mutating ``matplotlib.rcParams``).  If ``False``
        (default), the theme is applied only for the duration of each plot
        call, leaving global state untouched.
    **rc_overrides
        Individual ``rcParams`` entries to override in the theme.
    """

    _RASTER: frozenset[str] = frozenset({"png", "jpg", "jpeg", "tif", "tiff"})

    def __init__(
        self,
        palette: list[str] | None = None,
        global_theme: bool = False,
        **rc_overrides: Any,
    ) -> None:
        self._palette = palette if palette is not None else LOGO_PALETTE
        self._rc = theme_rc(self._palette, **rc_overrides)
        self._global = global_theme
        if global_theme:
            apply_global_theme(self._palette, **rc_overrides)

    # -- theming -----------------------------------------------------------

    def _theme_scope(self):
        """Context manager that scopes the theme unless it is applied globally."""
        if self._global:
            # Already applied globally; no scoping needed.
            from contextlib import nullcontext

            return nullcontext()
        return theme_context(self._rc)

    def _apply_plot_theme(self, plot: so.Plot) -> so.Plot:
        """Attach the house rc theme to a seaborn ``Plot`` (objects color cycle
        is independent of rcParams, so it is set separately via ``.scale``)."""
        return plot.theme(self._rc)

    # -- output ------------------------------------------------------------

    def __finish(
        self,
        fig: plt.Figure,
        show: bool,
        save: bool,
        filename: str,
        dpi: int = 400,
    ) -> None:
        """Lay out, optionally save (format from extension), optionally show."""
        fig.tight_layout()

        # Save before show: plt.show() may close/clear the figure.
        if save:
            ext = Path(filename).suffix.lower().lstrip(".")
            kw: dict[str, Any] = {"bbox_inches": "tight"}
            if ext in self._RASTER:
                kw["dpi"] = dpi
            fig.savefig(filename, **kw)

        if show:
            plt.show()

    @staticmethod
    def __scale_axes(plot: so.Plot, xscale: str, yscale: str) -> so.Plot:
        """Apply log scales to a ``Plot`` where requested (linear is the default)."""
        scales: dict[str, str] = {}
        if xscale == "log":
            scales["x"] = "log"
        if yscale == "log":
            scales["y"] = "log"
        return plot.scale(**scales) if scales else plot

    # -- generic single line plot (back-compat) ----------------------------

    def plot(
        self,
        x: list | float | np.ndarray,
        y: list | float | np.ndarray,
        fig: plt.Figure | None = None,
        ax: plt.Axes | None = None,
        xlabel: str = "",
        ylabel: str = "",
        xscale: str = "linear",
        yscale: str = "linear",
        title: str = "",
        label: str = "",
        grid: bool = True,
        show: bool = True,
        save: bool = False,
        filename: str = "plot.png",
        **line_kw: Any,
    ) -> tuple[plt.Figure, plt.Axes]:
        """Generic single line plot.

        Retained for backward compatibility; :func:`jaff.plotting.plot_rates`
        is the preferred entry point for one or many curves.

        Parameters
        ----------
        x, y
            Data to plot.
        fig, ax
            Existing figure/axes to draw onto.  Created if ``None``.
        label
            Legend entry; a legend is drawn when non-empty.
        **line_kw
            Forwarded to :class:`seaborn.objects.Line`.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        df = _frames.line_frame(x, y)

        with self._theme_scope():
            if fig is None or ax is None:
                fig, ax = plt.subplots()

            plot = so.Plot(df, x="x", y="y").add(
                so.Line(color=self._palette[0], **line_kw)
            )
            plot = self.__scale_axes(plot, xscale, yscale)
            plot = plot.label(x=xlabel, y=ylabel, title=title)
            self._apply_plot_theme(plot).on(ax).plot()

            ax.grid(grid)
            despine(ax)
            if label and ax.lines:
                # Objects marks carry no legend label for a single group; attach
                # the entry to the drawn line directly.
                ax.lines[-1].set_label(label)
                ax.legend()

            self.__finish(fig, show, save, filename)

        return fig, ax

    # -- generic multi-series renderer -------------------------------------

    def render_series(
        self,
        df: pd.DataFrame,
        *,
        x_col: str,
        y_col: str,
        group_col: str,
        xlabel: str,
        ylabel: str,
        log_x: bool,
        log_y: bool,
        dynamic_x: bool = False,
        trim: bool = False,
        shade: bool | float = False,
        bands_df: pd.DataFrame | None = None,
        layout: str = "overlay",
        title: str = "",
        fig: plt.Figure | None = None,
        ax: plt.Axes | None = None,
        grid: bool = True,
        show: bool = True,
        save: bool = False,
        filename: str = "plot.png",
    ) -> tuple[plt.Figure, Any]:
        """Render a tidy long frame of labelled curves.

        The single home for multi-curve rendering.  Both
        :func:`jaff.plotting.plot_rates` and :func:`jaff.plotting.plot_xsecs`
        build a frame and delegate here.

        Parameters
        ----------
        df : pandas.DataFrame
            Long frame with columns *x_col*, *y_col*, *group_col* (one group
            per curve).
        x_col, y_col, group_col : str
            Column names for the x, y, and grouping/label columns.
        xlabel, ylabel : str
            Axis labels.
        log_x, log_y : bool
            Log-scale the respective axis.
        dynamic_x : bool, optional
            If ``True``, drop the x-axis to linear when the positive data spans
            less than one decade (cross-section behaviour).  Default ``False``.
        trim : bool, optional
            Tighten the x-axis to the positive data span.  Default ``False``.
        shade : bool or float, optional
            Shade under each curve (``True`` = default alpha; float = alpha).
        bands_df : pandas.DataFrame or None, optional
            Band-averaged bars to overlay (columns ``lower``/``upper``/``xsec``
            already in plot units).
        layout : str, optional
            ``"overlay"`` (all curves on one axes) or ``"subplots"`` (one
            stacked panel per group).  Default ``"overlay"``.
        title, fig, ax, grid, show, save, filename
            Standard rendering controls.

        Returns
        -------
        tuple[Figure, Axes | numpy.ndarray]
            For ``"overlay"`` the second item is the single axes; for
            ``"subplots"`` it is the array of per-group axes.
        """
        if layout not in ("overlay", "subplots"):
            raise ValueError(f"layout must be 'overlay' or 'subplots', got {layout!r}")

        draw_kw = dict(
            x_col=x_col,
            y_col=y_col,
            group_col=group_col,
            xlabel=xlabel,
            ylabel=ylabel,
            log_x=log_x,
            log_y=log_y,
            dynamic_x=dynamic_x,
            trim=trim,
            shade=shade,
        )

        with self._theme_scope():
            if layout == "subplots":
                labels = list(dict.fromkeys(df[group_col]))
                n = len(labels)
                fig, axes = plt.subplots(
                    n, 1, sharex=True, figsize=(6.4, 2.6 * n), squeeze=False
                )
                axes = axes[:, 0]
                for i, (a, label) in enumerate(zip(axes, labels)):
                    sub = df[df[group_col] == label].reset_index(drop=True)
                    self.__draw_series(
                        a,
                        sub,
                        set_xlabel=(i == n - 1),  # only the bottom panel
                        title=label,
                        grid=grid,
                        **draw_kw,
                    )
                    self.__draw_bands(a, bands_df)
                if title:
                    fig.suptitle(title)
                self.__finish(fig, show, save, filename)
                return fig, axes

            if fig is None or ax is None:
                fig, ax = plt.subplots()
            self.__draw_series(ax, df, set_xlabel=True, title=title, grid=grid, **draw_kw)
            self.__draw_bands(ax, bands_df)
            self.__finish(fig, show, save, filename)
            return fig, ax

    def __draw_series(
        self,
        ax: plt.Axes,
        df: pd.DataFrame,
        *,
        x_col: str,
        y_col: str,
        group_col: str,
        xlabel: str,
        ylabel: str,
        log_x: bool,
        log_y: bool,
        dynamic_x: bool,
        trim: bool,
        shade: bool | float,
        grid: bool,
        set_xlabel: bool,
        title: str,
    ) -> None:
        """Draw the labelled curves in *df* onto a single axes."""
        labels = list(dict.fromkeys(df[group_col]))
        multi = len(labels) > 1

        # Span drives both dynamic-log and trim; compute once if either needs it.
        lo = hi = None
        if dynamic_x or trim:
            lo, hi = _xsec.positive_span(df, x_col, y_col)
            if lo is None:
                lo, hi = _xsec.finite_span(df, x_col)
        eff_log_x = _xsec.use_log_x(log_x, lo, hi) if dynamic_x else log_x

        plot = so.Plot(df, x=x_col, y=y_col)
        if multi:
            # Cycle the palette so more curves than colours still render.
            colors = [self._palette[i % len(self._palette)] for i in range(len(labels))]
            plot = plot.add(so.Line(), color=group_col).scale(
                color=so.Nominal(colors, order=labels)
            )
        else:
            plot = plot.add(so.Line(color=self._palette[0]))

        scales: dict[str, str] = {}
        if eff_log_x:
            scales["x"] = "log"
        if log_y:
            scales["y"] = "log"
        if scales:
            plot = plot.scale(**scales)

        plot = plot.label(
            x=xlabel if set_xlabel else "",
            y=ylabel,
            title=title,
        )
        self._apply_plot_theme(plot).on(ax).plot()

        if multi:
            self.__vary_linewidths(ax, len(labels))

        ax.grid(grid)
        despine(ax)
        if trim and lo is not None and hi is not None and hi > lo:
            ax.set_xlim(*_xsec.padded_limits(lo, hi, eff_log_x))

        if shade:
            self.__shade(ax, df, labels, multi, x_col, y_col, group_col, alpha=shade)
        if not set_xlabel:
            ax.set_xlabel("")

    def __vary_linewidths(self, ax: plt.Axes, n: int) -> None:
        """Give each overlaid curve a distinct width, thinnest drawn in front.

        Mirrors the natural line-width variation seen in seaborn ``relplot``
        overlays: the first curve is drawn thick, later curves get
        progressively thinner and a higher z-order so the thin line sits on top
        of the thick one where they overlap.

        Assigned by draw order: at this point ``ax.lines`` holds exactly the
        ``n`` series lines (in the ``Nominal`` group order), before any shade
        fills or band bars are added.
        """
        base = float(self._rc.get("lines.linewidth", 2.0))
        widths = np.linspace(base * 1.6, base * 0.75, n)
        lines = ax.lines
        for i in range(min(n, len(lines))):
            lines[i].set_linewidth(widths[i])
            lines[i].set_zorder(2.0 + i)  # later = thinner = higher zorder = front

    def __shade(
        self,
        ax: plt.Axes,
        df: pd.DataFrame,
        labels: list[str],
        multi: bool,
        x_col: str,
        y_col: str,
        group_col: str,
        *,
        alpha: bool | float,
    ) -> None:
        """Fill the area under each curve down to the axis bottom."""
        a = 0.18 if alpha is True else float(alpha)
        base = ax.get_ylim()[0]
        for i, label in enumerate(labels):
            sub = df[df[group_col] == label]
            color = self._palette[i % len(self._palette)] if multi else self._palette[0]
            ax.fill_between(
                sub[x_col].to_numpy(),
                base,
                sub[y_col].to_numpy(),
                color=color,
                alpha=a,
                linewidth=0,
            )

    def __draw_bands(self, ax: plt.Axes, band_df: pd.DataFrame | None) -> None:
        """Overlay band-averaged cross sections as bars, clipping open bands."""
        if band_df is None or band_df.empty:
            return
        _, xmax = ax.get_xlim()
        ybottom = ax.get_ylim()[0]
        lower = band_df["lower"].to_numpy(dtype=float)
        upper = band_df["upper"].to_numpy(dtype=float)
        # Clip an open (inf) or over-wide top edge to the visible axis range.
        upper = np.where(np.isfinite(upper), upper, xmax)
        upper = np.minimum(upper, xmax)
        width = np.clip(upper - lower, a_min=0.0, a_max=None)
        ax.bar(
            lower,
            band_df["xsec"].to_numpy(dtype=float),
            width=width,
            bottom=ybottom,
            align="edge",
            facecolor="none",
            edgecolor="#555555",
            linewidth=1.1,
            alpha=0.9,
            zorder=1.5,
            label="Band average",
        )

    # -- cross-section plot (back-compat single-xsecs entry) ---------------

    def plot_xsec(
        self,
        xsecs: XsecsProps,
        processes: list[str] | None = None,
        layout: str = "overlay",
        fig: plt.Figure | None = None,
        ax: plt.Axes | None = None,
        energy_unit: str = "eV",
        xsec_unit: str = "cm^2",
        energy_log: bool = True,
        xsec_log: bool = True,
        trim: bool = True,
        shade: bool | float = False,
        show_bands: bool = False,
        bands: pd.DataFrame | None = None,
        title: str = "",
        grid: bool = True,
        show: bool = True,
        save: bool = False,
        filename: str = "xsec.png",
    ) -> tuple[plt.Figure, Any]:
        """Plot the cross sections of a single ``xsecs`` mapping on log-log axes.

        Retained for backward compatibility; :func:`jaff.plotting.plot_xsecs`
        is the preferred entry point (it accepts one or many reactions).

        Parameters
        ----------
        xsecs
            Mapping as returned by :func:`jaff.physics.get_xsec` -- carries
            ``photon_energy`` (eV) plus any of ``photo_absorption``,
            ``photodecay`` (all in cm^2).
        processes
            Subset of process keys to draw.  Default: every process with data.
        layout, energy_unit, xsec_unit, energy_log, xsec_log, trim, shade,
        show_bands, bands, title, grid, show, save, filename
            See :meth:`render_series` and :func:`jaff.plotting.plot_xsecs`.

        Returns
        -------
        tuple[Figure, Axes | numpy.ndarray]
        """
        xsecs = _sampled_xsecs(xsecs)
        energy = xsecs["photon_energy"]
        if energy is None:
            raise ValueError("xsecs has no 'photon_energy' data to plot.")

        if processes is None:
            processes = [k for k in _units.PROCESS_LABELS if xsecs.get(k) is not None]
        processes = [p for p in processes if xsecs.get(p) is not None]
        if not processes:
            raise ValueError("xsecs has no cross-section data to plot.")

        x = np.asarray(_units.convert_energy(np.asarray(energy), "eV", energy_unit))
        series = [
            (
                x,
                np.asarray(_units.convert_xsec(np.asarray(xsecs[k]), "cm2", xsec_unit)),
                _units.PROCESS_LABELS.get(k, k),
            )
            for k in processes
        ]
        df = _frames.long_frame(series, "energy", "xsec", "process")
        band_df = (
            _frames.band_frame(bands, energy_unit, xsec_unit)
            if show_bands and bands is not None
            else None
        )

        return self.render_series(
            df,
            x_col="energy",
            y_col="xsec",
            group_col="process",
            xlabel=_units.energy_label(energy_unit),
            ylabel=_units.xsec_label(xsec_unit),
            log_x=energy_log,
            log_y=xsec_log,
            dynamic_x=True,
            trim=trim,
            shade=shade,
            bands_df=band_df,
            layout=layout,
            title=title,
            fig=fig,
            ax=ax,
            grid=grid,
            show=show,
            save=save,
            filename=filename,
        )
