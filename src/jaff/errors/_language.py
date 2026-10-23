# ABOUTME: Exception raised for an unsupported code-generation language or format
# ABOUTME: Raised by jaff.codegen.Language and the Codegen methods that resolve one


class InvalidLanguageError(Exception):
    """
    Raised when a code-generation language or language format is not supported.

    :class:`~jaff.codegen.Language` raises it when ``Language(name)`` is given a
    name (or alias) that no registered language answers to, and when a language
    view is requested with an unknown bracket or matrix format. The
    :class:`~jaff.codegen.Codegen` methods that take a ``lang`` argument
    propagate it unchanged.
    """

    pass
