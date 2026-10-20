"""
JAFF logging and progress-bar infrastructure.

This module sets up a single :class:`JaffLogger` singleton backed by Rich's
:class:`~rich.logging.RichHandler` (or a plain :class:`~logging.StreamHandler`
inside Jupyter notebooks).  It also provides :class:`JaffProgress`, a subclass
of :class:`~rich.progress.Progress` that adapts its rendering strategy based on
whether the process is running inside a Jupyter kernel, another interactive
Python session, or a plain terminal.

Module-level singletons
-----------------------
``jaff_console``
    The shared :class:`~rich.console.Console` instance used for all Rich output.
``jaff_progress``
    The global :class:`JaffProgress` bar.  It is started automatically at
    import time and registered with :func:`atexit.register` to stop cleanly
    when the interpreter exits (except in Jupyter / interactive sessions where
    the lifecycle is managed differently).
"""

import atexit
import logging
import os
import sys
from contextlib import contextmanager
from typing import Iterable, Optional, Sequence

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    ProgressType,
    SpinnerColumn,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
)
from rich.text import Text


def _is_jupyter() -> bool:
    """
    Return ``True`` when running inside a Jupyter (ZMQ) kernel.

    Detection is based on the class name of the active IPython shell; a
    :exc:`NameError` (``get_ipython`` not defined) means we are in a plain
    Python interpreter.

    Returns
    -------
    bool
        ``True`` iff the current process is a Jupyter kernel.
    """
    try:
        shell = get_ipython().__class__.__name__  # type: ignore
        return shell == "ZMQInteractiveShell"
    except NameError:
        return False


def _is_marimo() -> bool:
    """
    Return ``True`` when running inside a marimo notebook UI.

    Detection uses :func:`marimo.running_in_notebook`, which is ``True`` only
    when executing under the marimo editor/server (``marimo edit`` or
    ``marimo run``).  Running an exported marimo notebook as a plain script
    (``python notebook.py`` -> ``app.run()``) returns ``False``, so progress
    rendering falls back to the terminal strategy in that case.

    Returns
    -------
    bool
        ``True`` iff the current process is rendering inside a marimo notebook.
    """
    try:
        import marimo as mo

        return mo.running_in_notebook()
    except ImportError:
        return False


IN_JUPYTER = _is_jupyter()
IN_MARIMO = _is_marimo()
_IN_JUPYTER_ENV = IN_JUPYTER or "JPY_PARENT_PID" in os.environ
_IN_INTERACTIVE = hasattr(sys, "ps1")

_PROGRESS_KWARGS = dict(
    expand=True,
    redirect_stdout=False,
    redirect_stderr=False,
    transient=True,
)


def _make_progress_columns():
    """
    Return the default tuple of Rich progress columns used across JAFF.

    Columns (in order): spinner, description text, bar, M-of-N counter,
    elapsed time.

    Returns
    -------
    tuple
        A tuple of :class:`~rich.progress.ProgressColumn` instances.
    """
    return (
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=None),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
    )


def _make_progress(console: Console, auto_refresh: bool = True) -> Progress:
    """
    Construct a standalone :class:`~rich.progress.Progress` bar.

    Used internally to create a *temporary* progress context (e.g. inside
    Jupyter / interactive sessions) separate from the global
    :data:`jaff_progress` singleton.

    Parameters
    ----------
    console : rich.console.Console
        The Rich console to render into.
    auto_refresh : bool, optional
        Whether the progress bar refreshes automatically (default ``True``).

    Returns
    -------
    rich.progress.Progress
        A new Progress instance configured with the default JAFF columns.
    """
    return Progress(
        *_make_progress_columns(),
        console=console,
        auto_refresh=auto_refresh,
        **_PROGRESS_KWARGS,
    )


class JaffProgress(Progress):
    """
    Environment-aware Rich progress bar for JAFF.

    Extends :class:`~rich.progress.Progress` with two improvements:

    * :meth:`indeterminate` -- a context manager for tasks of unknown length
      that gracefully degrades in Jupyter and CI environments.
    * :meth:`track` -- overrides the parent implementation to handle the
      Jupyter / interactive-shell rendering quirks (Jupyter does not support
      in-place terminal overwrites, so a fresh :class:`~rich.progress.Progress`
      context is spun up for each iteration instead of mutating the global bar).

    In non-interactive terminal sessions the global progress bar is kept alive
    for the lifetime of the process via :func:`atexit.register`.
    """

    @contextmanager
    def indeterminate(self, description: str):
        """
        Context manager for a spinner task with unknown total work.

        Adapts to the execution environment:

        * **marimo notebook** -- uses :func:`marimo.status.spinner` so the
          spinner renders natively in the marimo UI.
        * **Jupyter / interactive** -- creates a temporary :class:`Progress`
          context so each indeterminate task gets its own isolated bar.
        * **Jupyter env without a kernel** (e.g. CI with ``JPY_PARENT_PID``) --
          yields ``None`` and does nothing.
        * **Terminal** -- adds the task to the shared global progress bar and
          removes it on exit.

        Parameters
        ----------
        description : str
            Human-readable label shown next to the spinner.

        Yields
        ------
        TaskID or None
            The Rich task identifier, or ``None`` in silent environments.
        """
        if IN_MARIMO:
            import marimo as mo

            with mo.status.spinner(title=description):
                yield None
            return

        if IN_JUPYTER or _IN_INTERACTIVE:
            with _make_progress(self.console) as p:
                task_id = p.add_task(description, total=None)
                yield task_id
        elif _IN_JUPYTER_ENV:
            yield None
        else:
            task_id = self.add_task(description, total=None)
            try:
                yield task_id
            finally:
                self.remove_task(task_id)

    def track(
        self,
        sequence: Sequence[ProgressType] | Iterable[ProgressType],
        total: Optional[float] = None,
        task_id: Optional[TaskID] = None,
        *args,
        description: str = "Working...",
        update_period: float = 0.1,
        **kwargs,
    ) -> Iterable[ProgressType]:
        """
        Iterate over *sequence* while rendering a progress bar.

        Overrides :meth:`rich.progress.Progress.track` to handle several
        environments:

        * **marimo notebook** -- delegates to :func:`marimo.status.progress_bar`
          so the bar renders natively in the marimo UI.
        * **Jupyter / interactive** -- creates a temporary
          :class:`~rich.progress.Progress` with ``auto_refresh=False`` and
          throttles refreshes to *update_period* seconds.
        * **Jupyter env (no kernel)** -- silently yields each value with no
          rendering.
        * **Terminal** -- uses the global shared progress bar, cleaning up
          finished tasks before adding a new one.

        Parameters
        ----------
        sequence : Sequence or Iterable
            The collection to iterate over.
        total : float or None, optional
            Total number of items.  Inferred from ``len(sequence)`` when
            ``None`` and the sequence supports it.
        task_id : TaskID or None, optional
            Existing task ID to update instead of creating a new one.
        description : str, optional
            Label displayed next to the progress bar (default ``"Working..."``).
        update_period : float, optional
            Minimum number of seconds between display refreshes
            (default ``0.1``).
        **kwargs
            Additional keyword arguments (accepted but unused; present for
            API compatibility with the parent class).

        Yields
        ------
        ProgressType
            Each item from *sequence* in order.
        """
        if IN_MARIMO:
            import marimo as mo

            if total is None and hasattr(sequence, "__len__"):
                total = len(sequence)  # type: ignore
            if total is not None:
                yield from mo.status.progress_bar(
                    sequence, title=description, total=total
                )
            else:
                with mo.status.spinner(title=description):
                    yield from sequence
            return

        if IN_JUPYTER or _IN_INTERACTIVE:
            with _make_progress(self.console, auto_refresh=False) as p:
                if total is None and hasattr(sequence, "__len__"):
                    total = len(sequence)  # type: ignore
                task_id = p.add_task(description, total=total)
                _get_time = p.get_time
                last_time = _get_time()
                for value in sequence:
                    yield value
                    p.advance(task_id, 1)
                    current_time = _get_time()
                    if (current_time - last_time) > update_period:
                        p.refresh()
                        last_time = current_time
                p.update(task_id, completed=total)
                p.refresh()
            return

        if _IN_JUPYTER_ENV:
            for value in sequence:
                yield value
            return

        for task in list(self.tasks):
            if task.finished:
                self.remove_task(task.id)

        if total is None and hasattr(sequence, "__len__"):
            total = len(sequence)  # type: ignore
        task_id = self.add_task(description, total=total)
        _finished = False
        try:
            _get_time = self.get_time
            last_time = _get_time()
            for value in sequence:
                yield value
                current_time = _get_time()
                if (current_time - last_time) > update_period:
                    self.advance(task_id, 1)
                    self.refresh()
                    last_time = current_time
                else:
                    self.advance(task_id, 1)
            _finished = True
            self.update(task_id, completed=total)
            self.refresh()
        finally:
            if not _finished:
                self.remove_task(task_id)


jaff_console = Console(force_jupyter=IN_JUPYTER)

jaff_progress = JaffProgress(
    *_make_progress_columns(), console=jaff_console, **_PROGRESS_KWARGS
)

if not _IN_JUPYTER_ENV and not _IN_INTERACTIVE and not IN_MARIMO:
    jaff_progress.start()
    atexit.register(jaff_progress.stop)


class _StripMarkupFormatter(logging.Formatter):
    """
    Logging formatter that strips Rich markup tags from messages.

    Used as the formatter for the plain :class:`~logging.StreamHandler` that
    handles output in Jupyter environments, where Rich's terminal escape codes
    would appear as raw text.
    """

    def format(self, record: logging.LogRecord) -> str:
        """
        Format *record*, removing any Rich markup before serialization.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to format.

        Returns
        -------
        str
            The formatted log string with Rich markup stripped.
        """
        original_msg, original_args = record.msg, record.args
        record.msg = Text.from_markup(record.getMessage()).plain
        record.args = None
        result = super().format(record)
        record.msg, record.args = original_msg, original_args
        return result


class JaffLogger:
    """
    Factory and wrapper for the named JAFF :class:`~logging.Logger`.

    Configures the logger on first access with either a
    :class:`~rich.logging.RichHandler` (terminal) or a plain
    :class:`~logging.StreamHandler` with markup stripping (Jupyter).
    Subsequent instantiations with the same *name* reuse the existing handler
    set because :func:`logging.getLogger` returns the same object.

    Parameters
    ----------
    name : str, optional
        Logger name (default ``"JAFF"``).
    level : int, optional
        Initial log level (default :data:`logging.INFO`).
    """

    def __init__(self, name: str = "JAFF", level: int = logging.INFO):
        """Configure and cache the named JAFF logger.

        Parameters
        ----------
        name : str, optional
            Logger name, by default ``"JAFF"``.
        level : int, optional
            Initial log level, by default :data:`logging.INFO`.
        """
        self.logger = logging.getLogger(name)

        if not self.logger.handlers:
            if _IN_JUPYTER_ENV or IN_MARIMO:
                handler: logging.Handler = logging.StreamHandler(sys.stdout)
                handler.setFormatter(_StripMarkupFormatter("%(levelname)-8s %(message)s"))
            else:
                handler = RichHandler(
                    console=jaff_console,
                    rich_tracebacks=True,
                    markup=True,
                    show_time=False,
                )
                handler.setFormatter(logging.Formatter("%(message)s"))
            self.logger.addHandler(handler)

        self.logger.setLevel(level)
        self.logger.propagate = False

    def get_logger(self) -> logging.Logger:
        """
        Return the underlying :class:`~logging.Logger` instance.

        Returns
        -------
        logging.Logger
            The configured logger.
        """
        return self.logger


def _route_warnings_through_logger() -> None:
    """Redirect :mod:`warnings` output through the JAFF logging stack.

    By default the :mod:`warnings` module writes straight to ``sys.stderr``,
    which tears the Rich live progress display (e.g. SciPy's
    ``IntegrationWarning`` emitted during rate-table integration).  Calling
    :func:`logging.captureWarnings` reroutes every warning to the
    ``py.warnings`` logger; binding that logger to a handler on
    :data:`jaff_console` lets Rich render warnings cleanly above the progress
    bar instead of corrupting it.  Warnings are still shown -- only the output
    sink changes.

    Idempotent: a second call is a no-op once the handler is attached.
    """
    logging.captureWarnings(True)
    warn_logger = logging.getLogger("py.warnings")
    if warn_logger.handlers:
        return

    if _IN_JUPYTER_ENV or IN_MARIMO:
        handler: logging.Handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(_StripMarkupFormatter("%(levelname)-8s %(message)s"))
    else:
        # markup=False: warning text may contain "[...]" that is not Rich markup.
        handler = RichHandler(
            console=jaff_console,
            markup=False,
            show_time=False,
            rich_tracebacks=False,
        )
        handler.setFormatter(logging.Formatter("%(message)s"))

    warn_logger.addHandler(handler)
    warn_logger.setLevel(logging.WARNING)
    warn_logger.propagate = False


_route_warnings_through_logger()
