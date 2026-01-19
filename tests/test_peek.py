# Tests for the peek function
# Verifies BAM file metadata extraction including contigs and modifications

import json

import pynanalogue


class TestPeek:
    """Tests for pynanalogue.peek function"""

    def test_peek_simple_bam(self, simple_bam):
        """Test peek returns correct contigs and modifications for simple_bam"""
        result = pynanalogue.peek(str(simple_bam))

        # Verify result structure
        assert "contigs" in result
        assert "modifications" in result

        # Verify contigs
        expected_contigs = {"contig_00000": 10000, "contig_00001": 10000}
        assert result["contigs"] == expected_contigs

        # Verify modifications
        expected_modifications = [["T", "+", "T"]]
        assert result["modifications"] == expected_modifications

    def test_peek_two_mods(self, two_mods_bam):
        """Test peek with BAM containing two different modifications"""
        result = pynanalogue.peek(str(two_mods_bam))

        # Verify result structure
        assert "contigs" in result
        assert "modifications" in result

        # Verify contigs
        expected_contigs = {"contig_00000": 10000, "contig_00001": 10000}
        assert result["contigs"] == expected_contigs

        # Verify two modifications detected (order may vary)
        assert len(result["modifications"]) == 2
        mods_set = {tuple(m) for m in result["modifications"]}
        assert ("T", "-", "T") in mods_set
        assert ("C", "+", "76792") in mods_set

    def test_peek_no_mods(self, tmp_path):
        """Test peek with BAM containing no modifications"""
        config = {
            "contigs": {"number": 3, "len_range": [20000, 20000]},
            "reads": [
                {
                    "number": 100,
                    "mapq_range": [20, 30],
                    "base_qual_range": [20, 30],
                    "len_range": [0.5, 0.5],
                    "insert_middle": "ATCG",
                    "mods": [],
                }
            ],
        }

        bam_path = tmp_path / "no_mods.bam"
        fasta_path = tmp_path / "no_mods.fasta"
        pynanalogue.simulate_mod_bam(
            json_config=json.dumps(config),
            bam_path=str(bam_path),
            fasta_path=str(fasta_path),
        )

        result = pynanalogue.peek(str(bam_path))

        # Verify result structure
        assert "contigs" in result
        assert "modifications" in result

        # Verify contigs
        expected_contigs = {
            "contig_00000": 20000,
            "contig_00001": 20000,
            "contig_00002": 20000,
        }
        assert result["contigs"] == expected_contigs

        # Verify no modifications detected
        assert result["modifications"] == []
