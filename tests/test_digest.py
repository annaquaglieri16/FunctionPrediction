import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from plm_benchmark.digest import digest, digest_multi_enzyme


def test_trypsin_respects_proline_rule():
    # K at position 3 (1-based) is followed by P -> NOT cleaved.
    # R at position 5 is followed by A -> cleaved.
    # K at position 8 is followed by A -> cleaved.
    # R at position 11 is the C-terminus -> cleaved.
    seq = "AAKPRAAKAAR"
    peptides = digest(seq, enzyme="trypsin", missed_cleavages=0, min_length=1)
    seqs = [p.sequence for p in peptides]
    assert seqs == ["AAKPR", "AAK", "AAR"]


def test_trypsin_missed_cleavages():
    seq = "AAKPRAAKAAR"
    peptides_1mc = {p.sequence for p in digest(seq, "trypsin", missed_cleavages=1, min_length=1)}
    assert peptides_1mc == {"AAKPR", "AAK", "AAR", "AAKPRAAK", "AAKAAR"}

    peptides_2mc = {p.sequence for p in digest(seq, "trypsin", missed_cleavages=2, min_length=1)}
    assert seq in peptides_2mc  # full sequence recovered with enough missed cleavages


def test_length_filters():
    seq = "AAKPRAAKAAR"
    peptides = digest(seq, "trypsin", missed_cleavages=2, min_length=4, max_length=6)
    for p in peptides:
        assert 4 <= len(p) <= 6


def test_lysc_cleaves_after_k_only():
    seq = "AKAKAR"
    peptides = [p.sequence for p in digest(seq, "lysc", missed_cleavages=0, min_length=1)]
    assert peptides == ["AK", "AK", "AR"]


def test_lysc_is_not_blocked_by_proline():
    # Unlike trypsin, Lys-C cleaves K-P: the ExPASy rule for lysc is a bare
    # "K" with no proline exception. Pinned explicitly because the sequence in
    # test_lysc_cleaves_after_k_only contains no K-P and so cannot see this.
    seq = "AAKPAAKAAR"
    peptides = [p.sequence for p in digest(seq, "lysc", missed_cleavages=0, min_length=1)]
    assert peptides == ["AAK", "PAAK", "AAR"]

    # Contrast: trypsin *is* blocked by proline at the same position.
    tryptic = [p.sequence for p in digest(seq, "trypsin", missed_cleavages=0, min_length=1)]
    assert tryptic == ["AAKPAAK", "AAR"]


def test_argc_cleaves_after_r_only():
    seq = "ARARAK"
    peptides = [p.sequence for p in digest(seq, "arg-c", missed_cleavages=0, min_length=1)]
    assert peptides == ["AR", "AR", "AK"]


def test_argc_is_not_blocked_by_proline():
    # Like Lys-C and unlike trypsin, the ExPASy arg-c rule is a bare "R" with
    # no proline exception. Pinned separately because ARARAK has no R-P.
    seq = "AARPAARAAK"
    peptides = [p.sequence for p in digest(seq, "arg-c", missed_cleavages=0, min_length=1)]
    assert peptides == ["AAR", "PAAR", "AAK"]


def test_peptide_positions_reconstruct_original():
    seq = "MKPRAAKRVGGSEQKPEPTIDEKAA"
    peptides = digest(seq, "trypsin", missed_cleavages=0, min_length=1)
    for p in peptides:
        assert seq[p.start:p.end] == p.sequence


def test_enzymes_config_names_resolve():
    # configs/enzymes.yaml names must exist in CLEAVAGE_RULES. Guards against
    # the rule table being renamed (e.g. swapped to ExPASy spellings) while
    # the config still refers to the old keys -- previously that failed only
    # at benchmark runtime, not in the tests.
    import yaml

    from plm_benchmark.digest import CLEAVAGE_RULES

    config_path = Path(__file__).resolve().parents[1] / "configs" / "enzymes.yaml"
    config = yaml.safe_load(config_path.read_text())

    names = [config["default_enzyme"], *config["enzymes_to_compare"]]
    missing = [n for n in names if n not in CLEAVAGE_RULES]
    assert not missing, f"unknown enzymes in enzymes.yaml: {missing}"


def test_unknown_enzyme_raises():
    with pytest.raises(ValueError):
        digest("AAKAAR", enzyme="not-a-real-enzyme")


def test_digest_multi_enzyme_is_union():
    seq = "AAKPRAAKAAR"
    trypsin_only = {p.sequence for p in digest(seq, "trypsin", missed_cleavages=0, min_length=1)}
    lysc_only = {p.sequence for p in digest(seq, "lysc", missed_cleavages=0, min_length=1)}
    combined = {
        p.sequence
        for p in digest_multi_enzyme(seq, ["trypsin", "lysc"], missed_cleavages=0, min_length=1)
    }
    assert trypsin_only <= combined
    assert lysc_only <= combined
