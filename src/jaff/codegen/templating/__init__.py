# ABOUTME: The $$JAFF template engine: grammar, single-pass block engine and commands
# ABOUTME: TemplateParser is the entry point; jaff.codegen re-exports it

from .functions import TemplateFunctions, jaffgen_function
from .parser import TemplateParser

__all__ = ["TemplateFunctions", "TemplateParser", "jaffgen_function"]
