# ABOUTME: Shielding plugins bound to H2 photodissociation (H2._PHOTON__H.H):
# ABOUTME: db1996 and hg2015 self-shielding
from .db1996 import DB1996
from .hg2015 import HG2015

__all__ = ["DB1996", "HG2015"]
