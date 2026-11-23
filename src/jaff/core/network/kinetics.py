# ABOUTME: Kinetics (net.kinetics): the chemistry's rates of change, i.e. the reaction
# ABOUTME: rates projected through the stoichiometry (species ODEs, weighted sums)
"""The chemical rates of change of a :class:`~jaff.Network`, ``net.kinetics``."""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

from sympy import Basic, Expr, Float, S

from ..reaction.rate import check_form

if TYPE_CHECKING:
    from .network import Network


class Kinetics:
    """The network's chemical rates of change.

    Reached as ``net.kinetics``.  Each quantity combines the reaction rates
    ``F_r`` (:meth:`Reactions.rates <jaff.core.Reactions.rates>`) with the
    stoichiometry ``ν_ri`` (products minus reactants of reaction *r* in species
    *i*): the species ODEs :meth:`odes`, and the rate of any weighted sum of the
    densities :meth:`weighted_rate`, with :meth:`dntot_dt` and :meth:`drho_dt` as
    its two common cases.

    Every builder takes ``form``: ``"expanded"`` writes each rate coefficient out,
    ``"reduced"`` uses its symbol ``krate[r]`` (see
    :class:`~jaff.core.reaction.Rate`).  Nothing is cached here: rate coefficients
    may still change after loading, so callers cache what they build.

    Attributes
    ----------
    net : Network
        The network whose reactions and species the rates describe.

    Examples
    --------
    >>> from jaff import Network
    >>> net = Network("h_photoionization")  # H + γ -> H+ + e-, H+ + e- -> H
    >>> [s.name for s in net.species.core]
    ['H', 'H+', 'e-']
    >>> odes = net.kinetics.odes()
    >>> odes[1] == -odes[0]
    True
    >>> net.kinetics.odes("reduced")[0]
    -krate[0]*nden[0] + krate[1]*nden[1]*nden[2]
    """

    def __init__(self, net: Network) -> None:
        """Bind to *net*; every quantity is built on request.

        Parameters
        ----------
        net : Network
            The network whose reactions and species the rates describe.
        """
        self.net: Network = net

    def odes(self, form: str = "expanded") -> list[Basic]:
        """Right-hand sides of the species ODEs, one per core species.

        ``d nden[s]/dt`` is the sum of the rates of the reactions producing *s*
        minus those destroying it; a species appearing several times on one side
        counts once per occurrence:

        ```
        d nden[s] / dt = sum_{r: s in products(r)} F_r
                       - sum_{r: s in reactants(r)} F_r
        ```

        Pseudo-species such as ``_PHOTON`` are not integrated.

        Parameters
        ----------
        form : {"expanded", "reduced"}, optional
            Form of the rates, by default ``"expanded"``.

        Returns
        -------
        list of sympy.Basic
            The time derivative of each core species' number density
            [cm⁻³ s⁻¹], in ``species.index`` order.

        Raises
        ------
        ValueError
            If *form* is unknown.
        """
        reactions = self.net.reactions
        rates = reactions.rates(form)
        odes: list[Basic] = [Float(0.0) for _ in range(self.net.species.core.count)]
        for rate, reaction in zip(rates, reactions, strict=True):
            for reactant in reaction.reactants.core:
                odes[reactant.index] -= rate
            for product in reaction.products.core:
                odes[product.index] += rate

        return odes

    def weighted_rate(
        self, weights: Sequence[Expr | float], rtol: float = 0.0, form: str = "expanded"
    ) -> Expr:
        """Rate of ``Σ_i w_i n_i``: ``Σ_r (Σ_i w_i ν_ri) F_r``.

        ``ν = product_matrix − reactant_matrix`` (core species).  The rate
        expressions are substituted, so the result can be differentiated and shares
        CSE with the species rows.  Reactions whose coefficient is exactly zero are
        skipped.

        Parameters
        ----------
        weights : Sequence[Expr | float]
            One weight per species, in species-index order.
        rtol : float, optional
            Numeric weights only: a reaction coefficient ``c_r = Σ_i w_i ν_ri`` is
            treated as zero when ``|c_r| <= rtol · Σ_i |w_i ν_ri|``, i.e. when it is a
            floating-point rounding residue of terms that cancel exactly in
            principle.  Default ``0.0`` keeps every non-zero coefficient.
        form : {"expanded", "reduced"}, optional
            Form of the rates, by default ``"expanded"``.

        Returns
        -------
        Expr
            The weighted rate [weight units · cm⁻³ s⁻¹].

        Raises
        ------
        ValueError
            If the number of weights differs from the number of species, or
            *form* is unknown.
        """
        net = self.net
        weights = list(weights)
        if len(weights) != net.species.count:
            raise ValueError(
                f"weighted_rate needs one weight per species "
                f"({net.species.count}), got {len(weights)}"
            )
        check_form(form)

        nu = net.product_matrix - net.reactant_matrix
        total = S.Zero
        for row, rate in zip(nu, net.reactions.rates(form)):
            terms = [w * int(n) for w, n in zip(weights, row) if n]
            coeff = sum(terms)
            if rtol and abs(coeff) <= rtol * sum(abs(t) for t in terms):
                continue

            if coeff != 0:
                total += coeff * rate

        return total

    def dntot_dt(self, form: str = "expanded") -> Expr:
        """Chemical rate of the total number density, ``Σ_r Δn_r F_r`` [cm⁻³ s⁻¹].

        :meth:`weighted_rate` with every weight ``1``.

        Parameters
        ----------
        form : {"expanded", "reduced"}, optional
            Form of the rates, by default ``"expanded"``.

        Returns
        -------
        Expr

        Raises
        ------
        ValueError
            If *form* is unknown.
        """
        return self.weighted_rate([1] * self.net.species.count, form=form)

    def drho_dt(self, form: str = "expanded") -> Expr:
        """Chemical rate of the mass density, ``Σ_r Δm_r F_r`` [g cm⁻³ s⁻¹].

        :meth:`weighted_rate` with the species masses as weights.  Never forced to
        zero: every reaction contributes its actual mass change, so the per-cell
        rate is kept even for reactions that pass :meth:`Reaction.check_mass`.  Only
        floating-point rounding residues (``|Δm_r| <= 1e-12 · Σ_i |m_i ν_ri|``, far
        below one electron mass) are dropped, so the expression -- and the Jacobian
        sparsity -- does not depend on the last bits of the species masses.  Unset
        masses count as ``0``.

        Parameters
        ----------
        form : {"expanded", "reduced"}, optional
            Form of the rates, by default ``"expanded"``.

        Returns
        -------
        Expr

        Raises
        ------
        ValueError
            If *form* is unknown.
        """
        masses = [s.mass or 0.0 for s in self.net.species]

        return self.weighted_rate(masses, rtol=1e-12, form=form)
