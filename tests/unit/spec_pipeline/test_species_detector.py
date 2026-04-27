"""Exhaustive-branch unit tests for detect_species(adata).

Covers every branch in the species_detector.py contract (locked):
  - Heuristic priority: declaration → ensembl_prefix → gene_case → fallthrough
  - Declaration normalization (scientific, common, genome builds; str + bytes)
  - Ensembl prefix ordering (ENSMUSG wins over ENSG longest-prefix-first)
  - Gene-case ambiguity (Title-case → mouse/rat, confidence=low)
  - Majority threshold 70% over first 100 var_names
  - SpeciesResult serialization (dict via asdict)
  - Merge semantics (overwrites only 'species' key in nvwa_meta)

Public API only: detect_species + SpeciesResult dataclass.
"""

from __future__ import annotations

import pytest

from src.domain.resolver.species_detector import (
    SpeciesResult,
    detect_species,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _set_var_names(adata, names: list[str]) -> None:
    """Replace adata.var_names with the given list of gene names.

    AnnData requires var.index length == n_vars. If names is shorter, pad with
    unique filler names. If names is longer, truncate to n_vars.
    """
    import pandas as pd
    n = adata.n_vars
    if len(names) < n:
        pad = [f"__pad__{i}" for i in range(n - len(names))]
        names = list(names) + pad
    else:
        names = list(names[:n])
    adata.var.index = pd.Index(names, name="var_names")


def _make_synthetic_adata(n_obs: int, var_names: list[str]):
    """Create a minimal synthetic AnnData with given var_names.

    Used for tests that need a specific number of genes (not the PBMC subset size).
    var_names becomes adata.var_names (the index of adata.var).
    """
    import anndata as ad
    import numpy as np
    import pandas as pd
    n_vars = len(var_names)
    X = np.zeros((n_obs, n_vars), dtype="float32")
    var_df = pd.DataFrame(index=pd.Index(var_names, name="var_names"))
    obs_df = pd.DataFrame(index=pd.RangeIndex(n_obs))
    return ad.AnnData(X=X, obs=obs_df, var=var_df)


def _make_names(template: list[str], pad_to: int, pad_gene: str) -> list[str]:
    """Pad a list of gene names to a target length using pad_gene."""
    result = list(template)
    while len(result) < pad_to:
        result.append(f"{pad_gene}{len(result)}")
    return result[:pad_to]


# ---------------------------------------------------------------------------
# 1. Declaration branch — priority: declaration wins everything else
# ---------------------------------------------------------------------------

class TestDeclarationBranch:

    def test_common_name_human(self, adata):
        adata.uns["species"] = "human"
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["confidence"] == "high"
        assert result["source"] == "declared"

    def test_scientific_name_homo_sapiens(self, adata):
        adata.uns["species"] = "Homo sapiens"
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["source"] == "declared"

    def test_genome_build_hg38(self, adata):
        adata.uns["species"] = "hg38"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_genome_build_hg19(self, adata):
        adata.uns["species"] = "hg19"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_genome_build_grch38(self, adata):
        adata.uns["species"] = "GRCh38"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_genome_build_grch37(self, adata):
        adata.uns["species"] = "GRCh37"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_common_name_mouse(self, adata):
        adata.uns["species"] = "mouse"
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "mouse"
        assert result["confidence"] == "high"
        assert result["source"] == "declared"

    def test_scientific_name_mus_musculus(self, adata):
        adata.uns["species"] = "Mus musculus"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    def test_genome_build_mm10(self, adata):
        adata.uns["species"] = "mm10"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    def test_genome_build_mm39(self, adata):
        adata.uns["species"] = "mm39"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    def test_genome_build_grcm38(self, adata):
        adata.uns["species"] = "GRCm38"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    def test_genome_build_grcm39(self, adata):
        adata.uns["species"] = "GRCm39"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    def test_common_name_rat(self, adata):
        adata.uns["species"] = "rat"
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "rat"
        assert result["confidence"] == "high"
        assert result["source"] == "declared"

    def test_scientific_name_rattus_norvegicus(self, adata):
        adata.uns["species"] = "Rattus norvegicus"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "rat"

    def test_genome_build_rn6(self, adata):
        adata.uns["species"] = "rn6"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "rat"

    def test_genome_build_rn7(self, adata):
        adata.uns["species"] = "rn7"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "rat"

    def test_genome_build_mratbn72(self, adata):
        adata.uns["species"] = "mRatBN7.2"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "rat"

    # --- case-insensitivity ---

    def test_declaration_case_insensitive_upper(self, adata):
        adata.uns["species"] = "HUMAN"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_declaration_case_insensitive_mixed(self, adata):
        adata.uns["species"] = "HoMo SaPiEnS"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_declaration_leading_trailing_whitespace(self, adata):
        adata.uns["species"] = "  mouse  "
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    # --- bytes input ---

    def test_declaration_as_bytes_human(self, adata):
        adata.uns["species"] = b"human"
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["source"] == "declared"

    def test_declaration_as_bytes_homo_sapiens(self, adata):
        adata.uns["species"] = b"Homo sapiens"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_declaration_as_bytes_mouse(self, adata):
        adata.uns["species"] = b"mouse"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    def test_declaration_as_bytes_rat(self, adata):
        adata.uns["species"] = b"rat"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "rat"

    # --- alternate key: 'organism' ---

    def test_organism_key_human(self, adata):
        adata.uns["organism"] = "human"
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["source"] == "declared"

    def test_organism_key_mouse(self, adata):
        adata.uns["organism"] = "Mus musculus"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "mouse"

    def test_organism_key_bytes(self, adata):
        adata.uns["organism"] = b"rat"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "rat"

    # --- species key wins over organism key (first in _DECLARATION_KEYS) ---

    def test_species_key_takes_priority_over_organism_key(self, adata):
        adata.uns["species"] = "human"
        adata.uns["organism"] = "mouse"
        detect_species(adata)
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    # --- unrecognized declaration falls through ---

    def test_unrecognized_declaration_falls_through(self, adata):
        """Unrecognized species declaration should not stop heuristic chain."""
        adata.uns["species"] = "zebrafish"
        # Give var_names that would resolve via ensembl prefix
        names = _make_names(
            [f"ENSG{i:011d}" for i in range(80)], 100, "UNKNOWN"
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        # Should fall through to ensembl prefix and find human
        assert result["source"] != "declared"

    def test_unrecognized_declaration_integer_ignored(self, adata):
        """Non-string, non-bytes declaration value is silently skipped."""
        adata.uns["species"] = 9606  # integer NCBI taxon ID for human
        # Give var_names that force a total fallthrough
        _set_var_names(adata, ["gene1", "gene2", "gene3"])
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        # Integer ignored → falls through → ultimately source='unknown'
        assert result["source"] != "declared"


# ---------------------------------------------------------------------------
# 2. Ensembl prefix branch
# ---------------------------------------------------------------------------

class TestEnsemblPrefixBranch:

    def test_human_ensg_majority(self, adata):
        """80/100 ENSG* genes → human, high confidence, ensembl_prefix source."""
        names = _make_names(
            [f"ENSG{i:011d}" for i in range(80)], 100, "MT-CO"
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["confidence"] == "high"
        assert result["source"] == "ensembl_prefix"

    def test_mouse_ensmusg_majority(self, adata):
        """80/100 ENSMUSG* genes → mouse, high confidence."""
        names = _make_names(
            [f"ENSMUSG{i:011d}" for i in range(80)], 100, "Mt-Co"
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "mouse"
        assert result["confidence"] == "high"
        assert result["source"] == "ensembl_prefix"

    def test_rat_ensrnog_majority(self, adata):
        """80/100 ENSRNOG* genes → rat, high confidence."""
        names = _make_names(
            [f"ENSRNOG{i:011d}" for i in range(80)], 100, "Mt-Co"
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "rat"
        assert result["confidence"] == "high"
        assert result["source"] == "ensembl_prefix"

    def test_ensmusg_wins_over_ensg_longest_prefix_first(self, adata):
        """ENSMUSG* IDs must not be double-counted as ENSG*.

        Mouse Ensembl IDs start with ENSMUSG; since ENSMUSG is a longer prefix
        than ENSG, the longest-first sort must assign them to 'mouse', not 'human'.
        """
        # All 100 are ENSMUSG — should → mouse, NOT human
        names = [f"ENSMUSG{i:011d}" for i in range(100)]
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "mouse", (
            "ENSMUSG prefix must win over ENSG (longest-prefix-first)"
        )
        assert result["source"] == "ensembl_prefix"

    def test_mixed_ensembl_majority_human_wins(self, adata):
        """75 ENSG + 25 ENSMUSG → human wins (≥70%)."""
        names = (
            [f"ENSG{i:011d}" for i in range(75)]
            + [f"ENSMUSG{i:011d}" for i in range(25)]
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["source"] == "ensembl_prefix"

    def test_mixed_ensembl_no_majority_falls_through(self, adata):
        """50 ENSG + 50 ENSMUSG — no species reaches 70%, falls through."""
        names = (
            [f"ENSG{i:011d}" for i in range(50)]
            + [f"ENSMUSG{i:011d}" for i in range(50)]
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] != "ensembl_prefix"

    def test_exactly_70_percent_threshold_passes(self, adata):
        """Exactly 70/100 Ensembl genes should meet the ≥0.70 threshold."""
        names = _make_names(
            [f"ENSG{i:011d}" for i in range(70)], 100, "MT-CO"
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] == "ensembl_prefix"
        assert result["value"] == "human"

    def test_below_70_percent_threshold_fails(self, adata):
        """69/100 Ensembl genes — just below threshold, should NOT match."""
        names = _make_names(
            [f"ENSG{i:011d}" for i in range(69)], 100, "MT-CO"
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] != "ensembl_prefix"

    def test_only_first_100_var_names_sampled(self, adata):
        """Even if var_names has 500 genes, only first 100 are sampled.

        Put 80 ENSG genes in positions 0-79, then non-Ensembl genes at 80-499.
        The 80 ENSG in the sample (80/100 = 80%) should still resolve to human.
        """
        import pandas as pd

        human_names = [f"ENSG{i:011d}" for i in range(80)]
        filler = [f"FILLER{i}" for i in range(420)]  # positions 80-499
        all_names = human_names + filler
        # Pad the current var_names to have at least 500 genes
        # We do this by creating a new adata with 500 var_names
        adata.var.index = pd.Index(all_names[:adata.n_vars], name="var_names")
        # Use first 100 logic — we need at least 100 genes for a clean test
        if adata.n_vars >= 100:
            detect_species(adata)
            result = adata.uns["nvwa_meta"]["species"]
            assert result["value"] == "human"
            assert result["source"] == "ensembl_prefix"

    def test_declaration_takes_priority_over_ensembl(self, adata):
        """Declaration branch fires before Ensembl prefix branch."""
        adata.uns["species"] = "mouse"  # declare mouse
        # Set ENSG genes (would say human via Ensembl)
        names = _make_names(
            [f"ENSG{i:011d}" for i in range(80)], 100, "MT-CO"
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "mouse"
        assert result["source"] == "declared"


# ---------------------------------------------------------------------------
# 3. Gene-case branch
# ---------------------------------------------------------------------------

class TestGeneCaseBranch:

    def test_all_upper_human_high_conf(self, adata):
        """80/100 all-uppercase gene symbols → human, high confidence."""
        names = _make_names(
            ["GAPDH", "IL2RA", "TP53", "ACTB", "MT-CO1",
             "CD3E", "CD8A", "FOXP3", "IL6", "TNF"] * 8,
            100,
            "GENE",
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["confidence"] == "high"
        assert result["source"] == "gene_case"

    def test_title_case_mouse_rat_ambiguity(self, adata):
        """80/100 Title-case gene symbols → candidates=[mouse,rat], confidence=low."""
        names = _make_names(
            ["Gapdh", "Trp53", "Actb", "Cd3e", "Foxp3",
             "Il6", "Tnf", "Bcl2", "Myc", "Tgfb1"] * 8,
            100,
            "Gene",
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] is None
        assert result["confidence"] == "low"
        assert result["source"] == "gene_case"
        assert result["candidates"] == ["mouse", "rat"]

    def test_title_case_candidates_exactly(self, adata):
        """Title-case ambiguity candidates must be exactly ['mouse', 'rat']."""
        names = ["Gapdh"] * 100
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        candidates = result.get("candidates")
        assert isinstance(candidates, list)
        assert set(candidates) == {"mouse", "rat"}
        assert len(candidates) == 2

    def test_gene_case_below_threshold_falls_through(self, adata):
        """60/100 uppercase genes (below 70%) — should not match gene_case."""
        upper = ["GAPDH"] * 60
        lower = ["gene_x"] * 40
        names = upper + lower
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] != "gene_case"

    def test_ensembl_priority_over_gene_case(self, adata):
        """Ensembl prefix branch fires before gene-case branch."""
        names = (
            [f"ENSG{i:011d}" for i in range(80)]
            + ["GAPDH"] * 20  # uppercase but Ensembl takes priority
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] == "ensembl_prefix"

    def test_declaration_priority_over_gene_case(self, adata):
        """Declaration fires before gene-case."""
        adata.uns["species"] = "rat"
        names = ["GAPDH"] * 100  # all-upper → would say human
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "rat"
        assert result["source"] == "declared"

    def test_mixed_case_no_majority_falls_through(self, adata):
        """Equal mix of upper and title-case — neither reaches 70%."""
        names = (["GAPDH"] * 50) + (["Gapdh"] * 50)
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] not in ("gene_case", "declared", "ensembl_prefix")

    def test_digits_and_hyphens_ignored_in_upper_check(self, adata):
        """Gene names with digits/hyphens are still all-upper if letters are upper."""
        names = ["MT-CO1", "IL2RA", "CD3E"] * 33 + ["MT-CO1"]  # all uppercase
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"
        assert result["source"] == "gene_case"


# ---------------------------------------------------------------------------
# 4. Total fallthrough branch
# ---------------------------------------------------------------------------

class TestFallthrough:

    def test_no_var_names_falls_through(self):
        """Empty var_names (synthetic AnnData) — all heuristics skip → unknown."""
        a = _make_synthetic_adata(n_obs=5, var_names=[])
        detect_species(a)
        result = a.uns["nvwa_meta"]["species"]
        assert result["value"] is None
        assert result["source"] == "unknown"
        assert result["confidence"] == "low"

    def test_fallthrough_candidates_all_species(self):
        """Fallthrough result must have candidates == _ALL_SPECIES."""
        a = _make_synthetic_adata(n_obs=5, var_names=[])
        detect_species(a)
        result = a.uns["nvwa_meta"]["species"]
        assert result["candidates"] == ["human", "mouse", "rat"]

    def test_ambiguous_gene_names_total_fallthrough(self, adata):
        """Purely ambiguous gene names that match no heuristic → unknown."""
        names = ["gene1", "gene2", "gene3"] * 34
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] == "unknown"
        assert result["value"] is None
        assert result["candidates"] == ["human", "mouse", "rat"]

    def test_no_declaration_no_ensembl_no_case(self, adata):
        """No declaration + mixed-case genes + no Ensembl IDs → unknown."""
        names = (["GAPDH"] * 30) + (["Gapdh"] * 30) + (["gene_x"] * 40)
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] == "unknown"


# ---------------------------------------------------------------------------
# 5. SpeciesResult serialization
# ---------------------------------------------------------------------------

class TestSpeciesResultSerialization:

    def test_result_stored_as_dict_not_dataclass(self, adata):
        """detect_species must store a plain dict, not a SpeciesResult object."""
        adata.uns["species"] = "human"
        detect_species(adata)
        stored = adata.uns["nvwa_meta"]["species"]
        assert isinstance(stored, dict), "Species result must be serialized to dict via asdict()"

    def test_dict_has_all_four_fields(self, adata):
        """Serialized dict must have value, confidence, source, candidates."""
        adata.uns["species"] = "human"
        detect_species(adata)
        stored = adata.uns["nvwa_meta"]["species"]
        assert "value" in stored
        assert "confidence" in stored
        assert "source" in stored
        assert "candidates" in stored

    def test_declared_human_candidates_is_none(self, adata):
        """Declaration path sets candidates=None (no ambiguity)."""
        adata.uns["species"] = "human"
        detect_species(adata)
        stored = adata.uns["nvwa_meta"]["species"]
        assert stored["candidates"] is None

    def test_fallthrough_candidates_is_list(self):
        """Fallthrough path sets candidates to a list of all species."""
        a = _make_synthetic_adata(n_obs=5, var_names=[])
        detect_species(a)
        stored = a.uns["nvwa_meta"]["species"]
        assert isinstance(stored["candidates"], list)

    def test_species_result_dataclass_fields(self):
        """SpeciesResult dataclass has exactly the documented fields."""
        import dataclasses
        fields = {f.name for f in dataclasses.fields(SpeciesResult)}
        assert fields == {"value", "confidence", "source", "candidates"}

    def test_species_result_default_candidates_none(self):
        """SpeciesResult candidates field defaults to None."""
        r = SpeciesResult(value="human", confidence="high", source="declared")
        assert r.candidates is None

    def test_species_result_candidates_list(self):
        """SpeciesResult can hold a list of candidates."""
        r = SpeciesResult(value=None, confidence="low", source="gene_case", candidates=["mouse", "rat"])
        assert r.candidates == ["mouse", "rat"]


# ---------------------------------------------------------------------------
# 6. Merge semantics — overwrites only 'species', preserves other nvwa_meta keys
# ---------------------------------------------------------------------------

class TestMergeSemantics:

    def test_preserves_schema_key(self, adata):
        """detect_species must not delete 'schema' from nvwa_meta."""
        adata.uns["nvwa_meta"] = {"schema": "v2", "other": "preserved"}
        adata.uns["species"] = "human"
        detect_species(adata)
        meta = adata.uns["nvwa_meta"]
        assert meta.get("schema") == "v2"

    def test_preserves_condition_cols_key(self, adata):
        """detect_species must not delete 'condition_cols' from nvwa_meta."""
        adata.uns["nvwa_meta"] = {"condition_cols": ["batch", "treatment"]}
        adata.uns["species"] = "mouse"
        detect_species(adata)
        meta = adata.uns["nvwa_meta"]
        assert meta.get("condition_cols") == ["batch", "treatment"]

    def test_preserves_ambiguous_cols_key(self, adata):
        """detect_species must not delete 'ambiguous_cols' from nvwa_meta."""
        adata.uns["nvwa_meta"] = {"ambiguous_cols": ["sex"]}
        adata.uns["species"] = "rat"
        detect_species(adata)
        meta = adata.uns["nvwa_meta"]
        assert meta.get("ambiguous_cols") == ["sex"]

    def test_preserves_user_resolutions_key(self, adata):
        """detect_species must not delete 'user_resolutions' from nvwa_meta."""
        adata.uns["nvwa_meta"] = {"user_resolutions": {"sex": "condition"}}
        adata.uns["species"] = "human"
        detect_species(adata)
        meta = adata.uns["nvwa_meta"]
        assert meta.get("user_resolutions") == {"sex": "condition"}

    def test_overwrites_existing_species_key(self, adata):
        """A pre-existing nvwa_meta['species'] dict is overwritten."""
        adata.uns["nvwa_meta"] = {
            "species": {"value": "mouse", "confidence": "low", "source": "unknown", "candidates": None},
            "schema": "v1",
        }
        adata.uns["species"] = "human"
        detect_species(adata)
        meta = adata.uns["nvwa_meta"]
        assert meta["species"]["value"] == "human"
        assert meta["species"]["source"] == "declared"
        # Other keys still intact
        assert meta["schema"] == "v1"

    def test_preserves_multiple_keys_simultaneously(self, adata):
        """All existing nvwa_meta keys survive; only 'species' is replaced."""
        adata.uns["nvwa_meta"] = {
            "schema": "v2",
            "condition_cols": ["batch"],
            "ambiguous_cols": [],
            "user_resolutions": {},
        }
        adata.uns["species"] = "mouse"
        detect_species(adata)
        meta = adata.uns["nvwa_meta"]
        assert meta["schema"] == "v2"
        assert meta["condition_cols"] == ["batch"]
        assert meta["ambiguous_cols"] == []
        assert meta["user_resolutions"] == {}
        assert meta["species"]["value"] == "mouse"

    def test_nvwa_meta_created_if_absent(self, adata):
        """detect_species creates nvwa_meta if it doesn't exist yet."""
        adata.uns.pop("nvwa_meta", None)
        adata.uns["species"] = "human"
        detect_species(adata)
        assert "nvwa_meta" in adata.uns
        assert adata.uns["nvwa_meta"]["species"]["value"] == "human"

    def test_idempotent_double_call(self, adata):
        """Calling detect_species twice gives same result; doesn't corrupt nvwa_meta."""
        adata.uns["species"] = "human"
        detect_species(adata)
        first = dict(adata.uns["nvwa_meta"])
        detect_species(adata)
        second = dict(adata.uns["nvwa_meta"])
        assert first == second

    def test_non_dict_nvwa_meta_replaced_gracefully(self, adata):
        """If nvwa_meta is not a dict, detect_species replaces it with a fresh dict."""
        adata.uns["nvwa_meta"] = "corrupted_string_value"
        adata.uns["species"] = "human"
        detect_species(adata)
        meta = adata.uns["nvwa_meta"]
        assert isinstance(meta, dict)
        assert meta["species"]["value"] == "human"


# ---------------------------------------------------------------------------
# 7. Edge cases and boundary conditions
# ---------------------------------------------------------------------------

class TestEdgeCases:

    def test_var_names_over_100_only_first_100_matter(self, adata):
        """With 150+ var_names, only first 100 are used for heuristics.

        Put 70 ENSG genes in positions 0-69 (just at threshold),
        then 80 ENSMUSG genes in positions 70-149.
        Should resolve to human (70/100 ENSG in sample = 70%).
        Uses _set_var_names so the index length always matches adata.n_vars.
        """
        # Build 150 test names, then let _set_var_names pad the rest to n_vars
        names_150 = (
            [f"ENSG{i:011d}" for i in range(70)]
            + [f"ENSMUSG{i:011d}" for i in range(80)]
        )
        # Only meaningful if we have >=100 genes in the fixture
        if adata.n_vars >= 100:
            _set_var_names(adata, names_150)
            detect_species(adata)
            result = adata.uns["nvwa_meta"]["species"]
            # First 100 = 70 ENSG + 30 ENSMUSG → human at exactly 70%
            assert result["source"] == "ensembl_prefix"
            assert result["value"] == "human"

    def test_fewer_than_100_var_names_uses_all(self):
        """If only 50 var_names exist, use all 50 for majority calculation.

        Uses a synthetic AnnData so we control n_vars exactly.
        """
        names = [f"ENSG{i:011d}" for i in range(40)] + [f"FILLER{i}" for i in range(10)]
        a = _make_synthetic_adata(n_obs=5, var_names=names)
        detect_species(a)
        result = a.uns["nvwa_meta"]["species"]
        # 40/50 = 80% ENSG → human
        assert result["value"] == "human"
        assert result["source"] == "ensembl_prefix"

    def test_single_var_name_ensg(self):
        """Single ENSG gene name (100% >= 70%) → human."""
        a = _make_synthetic_adata(n_obs=3, var_names=["ENSG00000000001"])
        detect_species(a)
        result = a.uns["nvwa_meta"]["species"]
        assert result["value"] == "human"

    def test_empty_uns(self):
        """Empty adata.uns with empty var_names — no crash, falls through to unknown."""
        a = _make_synthetic_adata(n_obs=3, var_names=[])
        a.uns.clear()
        detect_species(a)
        result = a.uns["nvwa_meta"]["species"]
        assert result["source"] == "unknown"

    def test_declaration_none_value_skipped(self):
        """adata.uns['species'] = None should be skipped (not treated as declaration)."""
        a = _make_synthetic_adata(n_obs=3, var_names=[])
        a.uns["species"] = None
        detect_species(a)
        result = a.uns["nvwa_meta"]["species"]
        assert result["source"] == "unknown"

    def test_gene_names_with_only_digits_not_title_case(self, adata):
        """Pure numeric gene names (e.g. '1234') should not count as title-case."""
        names = ["1234", "5678", "9012"] * 34
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        # No letter → not title-case, not upper → fallthrough
        assert result["source"] == "unknown"

    def test_ensmusg_prefix_not_misidentified_as_ensg(self, adata):
        """Regression: ENSMUSG must not score as ENSG (prefix collision)."""
        # 80 ENSMUSG IDs — should be mouse, not human
        names = [f"ENSMUSG{i:011d}" for i in range(80)] + [f"filler{i}" for i in range(20)]
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "mouse"
        assert result["value"] != "human"

    def test_ensrnog_not_misidentified_as_ensg(self, adata):
        """Regression: ENSRNOG must not score as ENSG (prefix collision)."""
        names = [f"ENSRNOG{i:011d}" for i in range(80)] + [f"filler{i}" for i in range(20)]
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["value"] == "rat"
        assert result["value"] != "human"

    def test_all_three_heuristics_skip_reaches_unknown(self, adata):
        """With no declaration, no Ensembl IDs, no clear case majority → unknown."""
        # Random-looking gene names: some upper, some lower, some numeric
        names = (
            ["GENA"] * 30  # upper
            + ["Genb"] * 30  # title-case
            + ["genc"] * 40  # all lower
        )
        _set_var_names(adata, names)
        detect_species(adata)
        result = adata.uns["nvwa_meta"]["species"]
        assert result["source"] == "unknown"
        assert result["value"] is None
        assert result["candidates"] == ["human", "mouse", "rat"]
