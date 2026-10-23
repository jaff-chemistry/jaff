# ABOUTME: Exception raised when a self-shielding function is registered twice
# ABOUTME: Used by the shielding-function registry in jaff.physics.photo_reactions


class RegistrationError(Exception):
    """
    Raised when a self-shielding function is registered under a key already in use.

    Shielding functions are registered by the key ``(name, reaction)``; a second
    class registered with the same ``name`` and ``reaction`` raises this error
    instead of silently replacing the first.
    """

    pass
