import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from plm_benchmark.labels import (
    EC_CLASS_NAMES,
    accession_from_header,
    load_ec_label_sets,
    load_ec_labels,
)

EXAMPLE_TSV = (
    "Entry\tEC number\tSubcellular location [CC]\tGene Ontology IDs\n"
    "P00350\t1.1.1.44\t\tGO:0004616; GO:0005829\n"
    "P0AD86\t\t\t\n"  # non-enzyme: no EC number
    "P00363\t1.3.5.1; 1.3.5.4\t\tGO:0009055\n"  # multifunctional: keep first
)


def test_load_ec_labels_parses_top_level_class(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    tsv_path.write_text(EXAMPLE_TSV)

    labels = load_ec_labels(tsv_path)

    assert labels == {"P00350": "1", "P00363": "1"}


def test_load_ec_labels_drops_proteins_without_ec(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    tsv_path.write_text(EXAMPLE_TSV)

    labels = load_ec_labels(tsv_path)

    assert "P0AD86" not in labels


def test_load_ec_labels_keeps_primary_ec_for_multifunctional_enzymes(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    tsv_path.write_text(EXAMPLE_TSV)

    labels = load_ec_labels(tsv_path)

    # "1.3.5.1; 1.3.5.4" -> primary is "1.3.5.1" -> top-level class "1"
    assert labels["P00363"] == "1"


def test_load_ec_labels_full_level_keeps_complete_ec_number(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    tsv_path.write_text(EXAMPLE_TSV)

    labels = load_ec_labels(tsv_path, level="full")

    assert labels == {"P00350": "1.1.1.44", "P00363": "1.3.5.1"}


def test_load_ec_labels_invalid_level_raises(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    tsv_path.write_text(EXAMPLE_TSV)

    with pytest.raises(ValueError):
        load_ec_labels(tsv_path, level="bogus")


def test_load_ec_label_sets_keeps_every_ec_number(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    tsv_path.write_text(EXAMPLE_TSV)

    label_sets = load_ec_label_sets(tsv_path)

    assert label_sets == {
        "P00350": frozenset({"1.1.1.44"}),
        "P00363": frozenset({"1.3.5.1", "1.3.5.4"}),
    }
    assert "P0AD86" not in label_sets  # non-enzyme, no EC


def test_load_ec_label_sets_differs_from_primary_only_parsing(tmp_path):
    tsv_path = tmp_path / "labels.tsv"
    tsv_path.write_text(EXAMPLE_TSV)

    primary = load_ec_labels(tsv_path, level="full")
    label_sets = load_ec_label_sets(tsv_path)

    # the secondary activity of the multifunctional enzyme is only visible in
    # the set version -- this is the label error the set form exists to fix
    assert "1.3.5.4" not in set(primary.values())
    assert "1.3.5.4" in label_sets["P00363"]


def test_accession_from_header_uniprot_style():
    assert accession_from_header("sp|P0AD86|LPT_ECOLI") == "P0AD86"
    assert accession_from_header("tr|A0A123|A0A123_ECOLI") == "A0A123"


def test_accession_from_header_passthrough_when_not_pipe_delimited():
    assert accession_from_header("P0AD86") == "P0AD86"


def test_ec_class_names_cover_all_seven_classes():
    assert set(EC_CLASS_NAMES) == {"1", "2", "3", "4", "5", "6", "7"}
