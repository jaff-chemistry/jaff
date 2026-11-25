# ABOUTME: jaff's built-in, opt-in plugins, imported as jaff.plugins: per-reaction models
# ABOUTME: (reactions/) and global physics (physics/radiation, physics/thermodynamics)
"""jaff plugins: self-contained, opt-in physics models, each selected by name in
``jaff.toml`` / ``jaffgen.toml``.

- ``reactions/<serialized reaction>/``: models bound to one reaction, by kind
  (``shielding/``, ...).
- ``physics/radiation/``: global radiation models (``shielding/``, ...).
- ``physics/thermodynamics/``: global heating / cooling mechanisms (their terms add up).

Every public subpackage is imported once, on the first plugin lookup, and its
``__init__.py`` imports the modules that register the plugins.  A ``.jfunc`` function
takes preference over a plugin of the same kind for a reaction.
"""
