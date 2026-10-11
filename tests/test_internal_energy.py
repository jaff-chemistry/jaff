# ABOUTME: Tests for InternalEnergy/DEDt arithmetic: same-class add/sub keeps the
# ABOUTME: class; mixed classes, plain operands and foreign networks are rejected

import operator

import pytest
from sympy import Symbol

from jaff.physics.thermodynamics.internal_energy import DEDt, InternalEnergy

TWO_SPECIES = "@format:idx,R,R,P,rate\n1,H,H,H2,1\n2,H2,H,H,1\n"
A, B = Symbol("a"), Symbol("b")


@pytest.fixture
def net(make_network):
    return make_network(TWO_SPECIES, funcfile=False)


@pytest.mark.parametrize("cls", [InternalEnergy, DEDt])
def test_add_same_class_keeps_class(net, cls):
    total = cls(A, net) + cls(B, net)
    assert type(total) is cls
    assert total.volumetric == A + B
    assert total._net is net


@pytest.mark.parametrize("cls", [InternalEnergy, DEDt])
def test_sub_same_class_keeps_class(net, cls):
    diff = cls(A, net) - cls(B, net)
    assert type(diff) is cls
    assert diff.volumetric == A - B


@pytest.mark.parametrize("op", [operator.add, operator.sub])
def test_mixed_classes_raise_type_error(net, op):
    energy, rate = InternalEnergy(A, net), DEDt(B, net)
    with pytest.raises(TypeError, match="'InternalEnergy' and 'DEDt'"):
        op(energy, rate)
    with pytest.raises(TypeError, match="'DEDt' and 'InternalEnergy'"):
        op(rate, energy)


@pytest.mark.parametrize("other", [1, 0.5, Symbol("x")])
@pytest.mark.parametrize("op", [operator.add, operator.sub])
def test_plain_operands_raise_type_error(net, other, op):
    energy = InternalEnergy(A, net)
    with pytest.raises(TypeError):
        op(energy, other)
    with pytest.raises(TypeError):
        op(other, energy)


def test_different_networks_raise_value_error(make_network):
    net_a = make_network(TWO_SPECIES, funcfile=False)
    net_b = make_network(TWO_SPECIES, funcfile=False)
    with pytest.raises(ValueError, match="different networks"):
        DEDt(A, net_a) + DEDt(B, net_b)
