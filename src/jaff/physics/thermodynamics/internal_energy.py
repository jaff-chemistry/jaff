# ABOUTME: InternalEnergy (symbolic internal energy in volumetric/specific/per-particle/
# ABOUTME: molar forms) and DEDt (its heating-rate counterpart with same-class +/-)

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, Callable, ClassVar

from sympy import Expr

from ..constants import N_A

if TYPE_CHECKING:
    from ...core import Network


class InternalEnergy:
    """Symbolic internal energy of a network in several normalisations.

    Wraps a volumetric internal energy ``E`` [erg cm⁻³] and derives the
    specific, per-particle and molar forms from the bound network's
    :attr:`~jaff.core.network._symbols.NetworkSymbols.rho` and
    :attr:`~jaff.core.network._symbols.NetworkSymbols.ntot`.
    """

    #: Form name -> accessor of the matching normalised property.
    _NORMALISED: ClassVar[dict[str, Callable[[InternalEnergy], Expr]]] = {
        "volumetric": lambda e: e.volumetric,
        "specific": lambda e: e.specific,
        "per_particle": lambda e: e.per_particle,
        "molar": lambda e: e.molar,
    }

    def __init__(self, expr: Expr, net: Network) -> None:
        """Wrap a volumetric internal energy.

        Parameters
        ----------
        expr : sympy.Expr
            Volumetric internal energy [erg cm⁻³].
        net : Network
            Network whose ``rho`` and ``ntot`` normalise *expr*.
        """
        self._net: Network = net
        self._vol_expr: Expr = expr

    def __add__(self, other: object) -> InternalEnergy:
        return self._combine(other, "+")

    def __sub__(self, other: object) -> InternalEnergy:
        return self._combine(other, "-")

    def __radd__(self, other: object) -> InternalEnergy:
        raise self._operand_error(other, "+", reflected=True)

    def __rsub__(self, other: object) -> InternalEnergy:
        raise self._operand_error(other, "-", reflected=True)

    def _combine(self, other: object, op: str) -> InternalEnergy:
        """Add or subtract another energy of exactly the same class and network.

        Parameters
        ----------
        other : object
            Right-hand operand.
        op : str
            ``"+"`` or ``"-"``.

        Returns
        -------
        InternalEnergy
            A new instance of ``type(self)`` bound to the same network.

        Raises
        ------
        TypeError
            If *other* is not exactly ``type(self)`` (subclasses do not mix).
        ValueError
            If *other* is bound to a different network.
        """
        if type(other) is not type(self):
            raise self._operand_error(other, op)

        assert isinstance(other, InternalEnergy)
        if other._net is not self._net:
            raise ValueError(f"Cannot {op} energies bound to different networks")

        if op == "+":
            expr = self._vol_expr + other._vol_expr
        else:
            expr = self._vol_expr - other._vol_expr

        return type(self)(expr, self._net)

    def _operand_error(
        self, other: object, op: str, reflected: bool = False
    ) -> TypeError:
        """``TypeError`` naming both operand types in evaluation order."""
        names = [type(self).__name__, type(other).__name__]
        left, right = reversed(names) if reflected else names

        return TypeError(f"unsupported operand type(s) for {op}: '{left}' and '{right}'")

    def normaliser(self, form: str) -> Expr:
        """Internal energy in the requested normalisation.

        Parameters
        ----------
        form : str
            ``"volumetric"``, ``"specific"``, ``"per_particle"`` or ``"molar"``.

        Returns
        -------
        sympy.Expr
            The internal energy normalised to *form* (see the property of the
            same name).

        Raises
        ------
        ValueError
            If *form* is not one of the forms above.
        """
        if form not in self._NORMALISED:
            raise ValueError(
                f"Invalid eos form: '{form}'. "
                f"Valid forms are: {', '.join(self._NORMALISED)}"
            )

        return self._NORMALISED[form](self)

    @cached_property
    def volumetric(self) -> Expr:
        """Internal energy per unit volume [erg cm⁻³]."""
        return self._vol_expr

    @cached_property
    def specific(self) -> Expr:
        """Internal energy per unit mass, ``E / ρ`` [erg g⁻¹]."""
        return self._vol_expr / self._net.symbols.rho

    @cached_property
    def per_particle(self) -> Expr:
        """Internal energy per particle, ``E / n_tot`` [erg]."""
        return self._vol_expr / self._net.symbols.ntot

    @cached_property
    def molar(self) -> Expr:
        """Internal energy per mole, ``N_A · E / n_tot`` [erg mol⁻¹]."""
        return self._vol_expr * N_A.cgs.value / self._net.symbols.ntot


class DEDt(InternalEnergy):
    """Rate of change of the volumetric internal energy, ``Ė`` [erg cm⁻³ s⁻¹].

    The normalised forms are time derivatives of the normalised energies, so they use
    the quotient rule with the network's per-cell density rates
    (:meth:`NetworkSymbols.drho_dt`, :meth:`NetworkSymbols.dntot_dt`, in this rate's
    :attr:`form`): neither ρ nor n_tot is constant in a cell.  ``E`` is read from
    ``net.thermodynamics.eos`` only when a normalised form is first requested, so the
    volumetric form never builds the EOS.

    Attributes
    ----------
    form : str
        ``"expanded"`` or ``"reduced"``: how the reaction rates in it are written
        (see :class:`~jaff.core.reaction.Rate`).  Only rates of the same form add.
    """

    def __init__(self, expr: Expr, net: Network, form: str = "expanded") -> None:
        """Wrap a volumetric energy rate.

        Parameters
        ----------
        expr : sympy.Expr
            Volumetric energy rate [erg cm⁻³ s⁻¹].
        net : Network
            Network whose ``rho`` and ``ntot`` normalise *expr*.
        form : {"expanded", "reduced"}, optional
            How the reaction rates in *expr* are written, by default
            ``"expanded"``.

        Raises
        ------
        ValueError
            If *form* is unknown.
        """
        from ...core.reaction.rate import check_form  # core imports physics

        check_form(form)
        super().__init__(expr, net)
        self.form: str = form

    def _combine(self, other: object, op: str) -> InternalEnergy:
        """Add or subtract another rate of the same network and form."""
        if isinstance(other, DEDt) and other.form != self.form:
            raise ValueError(f"Cannot {op} a {self.form} and a {other.form} energy rate")

        result = super()._combine(other, op)
        assert isinstance(result, DEDt)
        result.form = self.form

        return result

    @cached_property
    def _energy(self) -> Expr:
        """Volumetric internal energy ``E`` of the network's EOS."""
        return self._net.thermodynamics.eos.volumetric

    @cached_property
    def specific(self) -> Expr:
        """``d(E/ρ)/dt = Ė/ρ − E·ρ̇/ρ²`` [erg g⁻¹ s⁻¹]."""
        rho = self._net.symbols.rho

        return (
            self._vol_expr / rho
            - self._energy * self._net.symbols.drho_dt(self.form) / rho**2
        )

    @cached_property
    def per_particle(self) -> Expr:
        """``d(E/n)/dt = Ė/n − E·ṅ/n²`` [erg s⁻¹]."""
        ntot = self._net.symbols.ntot
        dntot = self._net.symbols.dntot_dt(self.form)

        return self._vol_expr / ntot - self._energy * dntot / ntot**2

    @cached_property
    def molar(self) -> Expr:
        """``N_A · d(E/n)/dt`` [erg mol⁻¹ s⁻¹]."""
        return self.per_particle * N_A.cgs.value
