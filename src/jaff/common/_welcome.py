"""
JAFF welcome banner (MOTD) generator.

Displays an ASCII-art title block with a randomly chosen adjective beginning
with the letter "F" -- because JAFF stands for *Just Another Fancy Format*
and the "F" word can be substituted for comedic effect.

Functions
---------
:func:`get_jaff_mode_dict`
    Return the title and sub-title strings for a given run mode.
:func:`motd`
    Return the full MOTD string for display at startup.
:func:`print_motd`
    Print the MOTD, degrading characters the console cannot encode.
"""

import sys
from textwrap import dedent

import numpy as np

words: list[str] = [
    "Fabulous",
    "Fair",
    "Faithful",
    "Familiar",
    "Famous",
    "Fancy",
    "Fantastic",
    "Fast",
    "Fearless",
    "Feisty",
    "Festive",
    "Fiery",
    "Fine",
    "Firm",
    "Fit",
    "Flaky",
    "Flamboyant",
    "Flashy",
    "Flat",
    "Flawed",
    "Flawless",
    "Flexible",
    "Flimsy",
    "Flippin",
    "Floccinaucinihilipilificating",
    "Flowery",
    "Fluffy",
    "Focused",
    "Formal",
    "Fortunate",
    "Fragile",
    "Friendly",
    "Frightening",
    "Frugal",
    "Functional",
    "Functionless",
    "Funky",
    "Fuzzy",
]


def get_jaff_mode_dict(fword: str, mode: str = "default") -> tuple[str, str]:
    """
    Return the ASCII-art title and sub-title strings for the requested run mode.

    Parameters
    ----------
    fword : str
        The "F word" to embed in the sub-title (e.g. ``"Fancy"``).  It will
        be title-cased before insertion.
    mode : {"default", "jaffgen", "jaffx"}, optional
        The run mode that determines which title block to use:

        * ``"default"`` -- plain ``JAFF`` banner.
        * ``"jaffgen"`` -- ``JAFFGEN`` code-generation banner.
        * ``"jaffx"`` -- ``JAFFX`` executor banner.

    Returns
    -------
    title : str
        Multi-line ASCII-art title string.
    sub_title : str
        One-line sub-title string, e.g. ``"Just Another Fancy Format!"``.
    """
    jaff_mode_dict: dict[str, dict[str, str]] = {
        "default": {
            "title": dedent(
                """
                ░░█ ▄▀█ █▀▀ █▀▀
                █▄█ █▀█ █▀░ █▀░
                """
            ),
            "sub_title": f"Just Another {fword.title()} Format!",
        },
        "jaffgen": {
            "title": dedent(
                """
                ░░█ ▄▀█ █▀▀ █▀▀ █▀▀ █▀▀ █▄░█
                █▄█ █▀█ █▀░ █▀░ █▄█ ██▄ █░▀█
                """
            ),
            "sub_title": f"Just Another {fword.title()} Format Generator!",
        },
        "jaffx": {
            "title": dedent(
                """
                ░░█ ▄▀█ █▀▀ █▀▀ ▀▄▀
                █▄█ █▀█ █▀░ █▀░ █░█
                """
            ),
            "sub_title": f"Just Another {fword.title()} Format Xecutor!",
        },
    }

    return jaff_mode_dict[mode]["title"], jaff_mode_dict[mode]["sub_title"]


def motd(mode: str = "default") -> str:
    """
    Generate the JAFF message-of-the-day banner string.

    Picks a random adjective from :data:`words` and assembles the title and
    sub-title using :func:`get_jaff_mode_dict`.

    Parameters
    ----------
    mode : {"default", "jaffgen", "jaffx"}, optional
        Run mode passed to :func:`get_jaff_mode_dict` (default ``"default"``).

    Returns
    -------
    str
        The complete MOTD string, ready to print to the terminal.
    """
    fword = np.random.choice(words)
    title, sub_title = get_jaff_mode_dict(fword, mode)

    return f"{title}\n{sub_title}\n"


def print_motd(mode: str = "default") -> None:
    """
    Print the JAFF banner to ``sys.stdout`` on any output encoding.

    The banner's block-art characters are not representable in every encoding:
    on Windows, a piped or redirected stdout uses the locale code page (e.g.
    cp1252), where printing them raises :class:`UnicodeEncodeError`.  Characters
    the stream cannot encode are replaced (``?``) instead.

    Parameters
    ----------
    mode : {"default", "jaffgen", "jaffx"}, optional
        Run mode passed to :func:`motd` (default ``"default"``).
    """
    text = motd(mode)
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"

    print(text.encode(encoding, errors="replace").decode(encoding))
