"""Thermodynamics backend abstraction.

The whole pipeline talks to RNA folding through a single interface so the
design logic never depends on which engine is installed.

* :class:`NupackBackend` mirrors the Toehold-VISTA reference model exactly
  (``Model(material='rna', ensemble='stacking', celsius=T, sodium, magnesium)``
  with ``mfe`` / ``defect`` / ``pfunc``).  Used automatically when the
  ``nupack`` module can be imported.
* :class:`ViennaRNABackend` reproduces the same quantities from ViennaRNA's
  partition function and base-pair-probability matrix.  This is the fallback
  and is what runs on a stock scientific-Python install.

Quantities exposed
------------------
mfe(seq)                     -> (dot_bracket, energy)
complex_mfe(seqs)            -> energy of the ordered complex (cofold)
ensemble_defect(seq, struct) -> normalised ensemble defect (0..1)
    * "SED" (specified) : pass an intended structure
    * "NED" (native)    : pass the MFE structure
accessibility(seq, i, L)     -> mean unpaired probability over [i, i+L)
unpaired_probabilities(seq)  -> per-position P(unpaired)
pair_probabilities(seq, thr) -> [(i, j, p), ...] 1-based, for arc plots
binding_dG(a, b)             -> Delta G of duplex formation (kcal/mol)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Sequence

from . import sequence_utils as su


# --------------------------------------------------------------------------- #
# structure helpers                                                           #
# --------------------------------------------------------------------------- #
def parse_pairs(structure: str) -> list[int]:
    """Return a 1-based partner array for a dot-bracket string.

    ``partner[i]`` is the 1-based index paired with i, or 0 if i is unpaired.
    Handles round/square/curly brackets and treats ``+``/``&`` strand breaks
    and other characters as unpaired placeholders.
    """
    n = len(structure)
    partner = [0] * (n + 1)
    stacks: dict[str, list[int]] = {"(": [], "[": [], "{": []}
    close_to_open = {")": "(", "]": "[", "}": "{"}
    for idx, ch in enumerate(structure, start=1):
        if ch in stacks:
            stacks[ch].append(idx)
        elif ch in close_to_open:
            op = close_to_open[ch]
            if stacks[op]:
                j = stacks[op].pop()
                partner[idx] = j
                partner[j] = idx
    return partner


class ThermoBackend:
    name = "abstract"

    def mfe(self, seq: str) -> tuple[str, float]:
        raise NotImplementedError

    def complex_mfe(self, seqs: Sequence[str]) -> float:
        raise NotImplementedError

    def ensemble_defect(self, seq: str, structure: str) -> float:
        raise NotImplementedError

    def unpaired_probabilities(self, seq: str) -> list[float]:
        raise NotImplementedError

    def pair_probabilities(self, seq: str, threshold: float = 0.01):
        raise NotImplementedError

    def ensemble_energy(self, seq: str) -> float:
        """Ensemble free energy -RT ln(Q) of the sequence's fold ensemble."""
        raise NotImplementedError

    def structure_probability(self, seq: str, temperature_c: float = 37.0) -> float:
        """Boltzmann probability that a sequence adopts its own MFE structure,
        p = exp(-(E_mfe - E_ensemble)/kT).  1.0 == the MFE fold dominates the
        ensemble; near 0 == the fold is unreliable."""
        _struct, e_mfe = self.mfe(seq)
        e_ens = self.ensemble_energy(seq)
        kT = 0.001987 * (temperature_c + 273.15)
        return math.exp(-(e_mfe - e_ens) / kT) if kT else 0.0

    def unpaired_probabilities_constrained(self, seq: str,
                                           forced_unpaired) -> list[float]:
        """Per-position P(unpaired) when the 0-based positions in
        ``forced_unpaired`` are held single-stranded (emulating a trigger that
        has hybridised there).  Backends without constraint support fall back to
        the unconstrained distribution."""
        return self.unpaired_probabilities(seq)

    # ---- quantities derived once the primitives above exist -------------- #
    def accessibility(self, seq: str, start: int, length: int) -> float:
        up = self.unpaired_probabilities(seq)
        window = up[start:start + length]
        return sum(window) / len(window) if window else 0.0

    def region_accessibility(self, seq: str, indices, forced_unpaired=None) -> float:
        """Mean unpaired probability over a set of 0-based positions, optionally
        conditioning on ``forced_unpaired`` positions being single-stranded."""
        if not indices:
            return 0.0
        up = (self.unpaired_probabilities_constrained(seq, forced_unpaired)
              if forced_unpaired else self.unpaired_probabilities(seq))
        vals = [up[i] for i in indices if 0 <= i < len(up)]
        return sum(vals) / len(vals) if vals else 0.0

    def native_defect(self, seq: str) -> float:
        struct, _ = self.mfe(seq)
        return self.ensemble_defect(seq, struct)

    def open_defect(self, seq: str) -> float:
        """SED against the fully single-stranded reference -- a direct
        'how structured is this region' score (0 == perfectly open)."""
        return self.ensemble_defect(seq, "." * len(seq))

    def binding_dG(self, a: str, b: str) -> float:
        """Delta G of forming the a:b duplex, ``G(ab) - G(a) - G(b)`` using
        complex/monomer MFE.  More negative == tighter binding."""
        g_ab = self.complex_mfe([a, b])
        g_a = self.mfe(a)[1]
        g_b = self.mfe(b)[1]
        return g_ab - g_a - g_b


# --------------------------------------------------------------------------- #
# ViennaRNA implementation                                                    #
# --------------------------------------------------------------------------- #
class ViennaRNABackend(ThermoBackend):
    name = "ViennaRNA"

    def __init__(self, temperature_c: float = 37.0):
        import RNA  # noqa: local import so the module loads without RNA present
        self._RNA = RNA
        self.T = float(temperature_c)
        RNA.cvar.temperature = self.T  # affects cofold() and legacy calls
        # dangles/params left at ViennaRNA defaults (Turner 2004), matching the
        # VISTA notebooks' ViennaRNA usage.

    def _md(self):
        md = self._RNA.md()
        md.temperature = self.T
        return md

    def _fc(self, seq: str):
        return self._RNA.fold_compound(su.to_rna(seq), self._md())

    def mfe(self, seq: str) -> tuple[str, float]:
        fc = self._fc(seq)
        ss, e = fc.mfe()
        return ss, float(e)

    def complex_mfe(self, seqs: Sequence[str]) -> float:
        # fold_compound + mfe() supports 2- and >=3-strand complexes (RNA.cofold
        # only handles two strands).
        joined = "&".join(su.to_rna(s) for s in seqs)
        fc = self._RNA.fold_compound(joined, self._md())
        _, e = fc.mfe()
        return float(e)

    def ensemble_energy(self, seq: str) -> float:
        fc = self._fc(seq)
        _struct, e = fc.pf()
        return float(e)

    def _bpp_and_unpaired(self, seq: str):
        rna = su.to_rna(seq)
        n = len(rna)
        fc = self._RNA.fold_compound(rna, self._md())
        fc.pf()
        bpp = fc.bpp()
        paired_sum = [0.0] * (n + 1)
        pairs: list[tuple[int, int, float]] = []
        for i in range(1, n + 1):
            row = bpp[i]
            for j in range(i + 1, n + 1):
                p = row[j]
                if p > 0.0:
                    paired_sum[i] += p
                    paired_sum[j] += p
                    pairs.append((i, j, p))
        unpaired = [max(0.0, 1.0 - paired_sum[i]) for i in range(1, n + 1)]
        return n, paired_sum, unpaired, pairs

    def ensemble_defect(self, seq: str, structure: str) -> float:
        rna = su.to_rna(seq)
        if len(structure) != len(rna):
            raise ValueError(
                f"structure length {len(structure)} != sequence length {len(rna)}")
        n, paired_sum, _unpaired, _pairs = self._bpp_and_unpaired(rna)
        # need the actual pair probability for specified pairs
        fc = self._RNA.fold_compound(rna, self._md())
        fc.pf()
        bpp = fc.bpp()
        partner = parse_pairs(structure)
        defect = 0.0
        for i in range(1, n + 1):
            j = partner[i]
            if j == 0:                       # specified unpaired
                p_correct = max(0.0, 1.0 - paired_sum[i])
            else:                            # specified paired i-j
                a, b = (i, j) if i < j else (j, i)
                p_correct = bpp[a][b]
            defect += (1.0 - p_correct)
        return defect / n

    def unpaired_probabilities(self, seq: str) -> list[float]:
        _n, _ps, unpaired, _pairs = self._bpp_and_unpaired(seq)
        return unpaired

    def unpaired_probabilities_constrained(self, seq: str,
                                           forced_unpaired) -> list[float]:
        rna = su.to_rna(seq)
        n = len(rna)
        fc = self._RNA.fold_compound(rna, self._md())
        for i in forced_unpaired:
            if 0 <= i < n:
                fc.hc_add_up(i + 1)   # ViennaRNA is 1-based
        fc.pf()
        bpp = fc.bpp()
        paired_sum = [0.0] * (n + 1)
        for i in range(1, n + 1):
            row = bpp[i]
            for j in range(i + 1, n + 1):
                p = row[j]
                if p > 0.0:
                    paired_sum[i] += p
                    paired_sum[j] += p
        return [max(0.0, 1.0 - paired_sum[i]) for i in range(1, n + 1)]

    def pair_probabilities(self, seq: str, threshold: float = 0.01):
        _n, _ps, _unpaired, pairs = self._bpp_and_unpaired(seq)
        return [(i, j, p) for (i, j, p) in pairs if p >= threshold]


# --------------------------------------------------------------------------- #
# NUPACK implementation (used when installed; mirrors Toehold-VISTA)          #
# --------------------------------------------------------------------------- #
class NupackBackend(ThermoBackend):
    name = "NUPACK"

    def __init__(self, temperature_c=37.0, sodium=1.0, magnesium=0.0,
                 material="rna", ensemble="stacking"):
        import nupack  # noqa
        self._nupack = nupack
        self._temp = float(temperature_c)
        self._vienna_helper = None
        self.model = nupack.Model(
            material=material, ensemble=ensemble,
            celsius=float(temperature_c), sodium=sodium, magnesium=magnesium)

    def mfe(self, seq: str) -> tuple[str, float]:
        res = self._nupack.mfe(strands=[su.to_rna(seq)], model=self.model)
        return str(res[0].structure), float(res[0].energy)

    def complex_mfe(self, seqs: Sequence[str]) -> float:
        res = self._nupack.mfe(strands=[su.to_rna(s) for s in seqs],
                               model=self.model)
        return float(res[0].energy)

    def ensemble_defect(self, seq: str, structure: str) -> float:
        # NUPACK 4's defect() already returns the NORMALISED ensemble defect
        # (fraction of incorrectly paired nucleotides, 0..1) -- do NOT divide by
        # length again (verified: nupack_defect * len ~= ViennaRNA normalised).
        return float(self._nupack.defect(strands=[su.to_rna(seq)],
                                         structure=structure, model=self.model))

    def ensemble_energy(self, seq: str) -> float:
        pf = self._nupack.pfunc(strands=[su.to_rna(seq)], model=self.model)
        return float(pf[1])

    def _pair_matrix(self, seq: str):
        rna = su.to_rna(seq)
        result = self._nupack.pairs(strands=[rna], model=self.model)
        return rna, result.to_array()

    def unpaired_probabilities(self, seq: str) -> list[float]:
        rna, mat = self._pair_matrix(seq)
        n = len(rna)
        # NUPACK 4.1 PairMatrix.to_array() is n x n with each row summing to 1;
        # the diagonal entry mat[i][i] is P(base i unpaired).
        return [float(mat[i][i]) for i in range(n)]

    def unpaired_probabilities_constrained(self, seq: str,
                                           forced_unpaired) -> list[float]:
        # NUPACK's analysis API has no hard-constraint hook for forcing bases
        # single-stranded; ViennaRNA (always installed) does. Unconstrained
        # accessibility agrees between the two engines to <0.01, so borrowing
        # ViennaRNA for this one quantity keeps the AND-mechanism scoring intact
        # without distorting results.
        if self._vienna_helper is None:
            self._vienna_helper = ViennaRNABackend(temperature_c=self._temp)
        return self._vienna_helper.unpaired_probabilities_constrained(
            seq, forced_unpaired)

    def pair_probabilities(self, seq: str, threshold: float = 0.01):
        rna, mat = self._pair_matrix(seq)
        n = len(rna)
        out = []
        for i in range(n):
            for j in range(i + 1, n):
                p = float(mat[i][j])
                if p >= threshold:
                    out.append((i + 1, j + 1, p))
        return out


# --------------------------------------------------------------------------- #
# factory                                                                     #
# --------------------------------------------------------------------------- #
@dataclass
class _BackendKey:
    prefer_nupack: bool
    temperature_c: float
    sodium: float
    magnesium: float
    material: str
    ensemble: str


def _make_backend(cfg) -> ThermoBackend:
    if getattr(cfg, "prefer_nupack", True):
        try:
            return NupackBackend(
                temperature_c=cfg.temperature_c, sodium=cfg.sodium,
                magnesium=cfg.magnesium, material=cfg.material,
                ensemble=cfg.ensemble)
        except Exception:
            pass
    return ViennaRNABackend(temperature_c=cfg.temperature_c)


_CACHE: dict = {}


def get_backend(cfg) -> ThermoBackend:
    """Return a cached backend for this config's physical model."""
    key = (bool(getattr(cfg, "prefer_nupack", True)), cfg.temperature_c,
           cfg.sodium, cfg.magnesium, cfg.material, cfg.ensemble)
    if key not in _CACHE:
        _CACHE[key] = _make_backend(cfg)
    return _CACHE[key]
