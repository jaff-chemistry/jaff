# ABOUTME: Registry of $$JAFF commands keyed by the grammar's command aliases
# ABOUTME: The syntax layer builds commands through it; one package per command

from typing import Dict, Type

from .base import Command, FileSetting
from .call import Call
from .comment import Comment
from .define import Define
from .end import End
from .get import Get
from .has import Has
from .lang import Lang
from .reduce import Reduce
from .repeat import Repeat
from .sub import Sub

COMMANDS: Dict[str, Type[Command]] = {
    "repeat": Repeat,
    "sub": Sub,
    "get": Get,
    "has": Has,
    "reduce": Reduce,
    "define": Define,
    "end": End,
    "lang": Lang,
    "comment": Comment,
    "call": Call,
}

__all__ = [
    "COMMANDS",
    "Call",
    "Command",
    "Comment",
    "Define",
    "End",
    "FileSetting",
    "Get",
    "Has",
    "Lang",
    "Reduce",
    "Repeat",
    "Sub",
]
