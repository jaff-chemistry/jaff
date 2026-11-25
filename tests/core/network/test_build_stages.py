# ABOUTME: Network build stages: species and elements are complete before the first
# ABOUTME: reaction is built, so reaction inputs can be resolved against all species

from jaff.core.reaction import Reaction

LINES = "H + H -> H2  []  1e-17\nH2 + C -> CH + H  []  1e-10\nCH -> C + H  []  1e-9\n"


def test_species_complete_before_first_reaction(make_network, monkeypatch):
    seen = []
    original = Reaction.__init__

    def spy(self, net, *args, **kwargs):
        seen.append(([s.name for s in net.species], net.elements))
        original(self, net, *args, **kwargs)

    monkeypatch.setattr(Reaction, "__init__", spy)
    net = make_network(LINES, name="s.jet")

    names, elements = seen[0]
    assert names == [s.name for s in net.species]
    assert elements is net.elements


def test_species_keep_file_order(make_network):
    net = make_network(LINES, name="s.jet")

    assert [s.name for s in net.species] == ["H", "H2", "C", "CH"]
