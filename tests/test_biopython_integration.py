import pytest
from unittest.mock import patch, mock_open
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from AlignmentViewer import AlignmentViewer
from AlignmentViewer.config import DisplayConfig


class TestBioPythonIntegration:
    """Integration tests for BioPython SeqRecord support"""

    def test_display_alignment_with_seqio_list(self):
        """Test the exact use case: sequences = list(SeqIO.parse(...)); AlignmentViewer.display_alignment(sequences, ...)"""
        # Create mock SeqRecord objects that would come from SeqIO.parse
        mock_seqrecords = [
            SeqRecord(Seq("ACGTACGTACGT"), id="seq1"),
            SeqRecord(Seq("TGCATGCATGCA"), id="seq2"),
            SeqRecord(Seq("GGCCGGCCGGCC"), id="seq3"),
            SeqRecord(Seq("AATTAATTAATT"), id="seq4"),
            SeqRecord(Seq("CCGGCCGGCCGG"), id="seq5"),
        ]

        # Test with as_html=True (the exact case from the notebook)
        html_result = AlignmentViewer.get_alignment_html(
            mock_seqrecords,
            nseqs=5,
            ncols=60,
            as_html=True
        )

        # Verify we get HTML output
        assert isinstance(html_result, str)
        assert "<div" in html_result
        assert "seq1" in html_result
        # Check that individual nucleotides are present (they're in separate spans)
        assert ">A</span>" in html_result  # First nucleotide from seq1
        assert ">C</span>" in html_result  # Second nucleotide from seq1

    def test_display_alignment_with_seqio_list_terminal(self):
        """Test with as_html=False"""
        mock_seqrecords = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("TGCA"), id="seq2"),
        ]

        # This should not raise an exception
        try:
            AlignmentViewer.display_alignment(
                mock_seqrecords,
                nseqs=2,
                ncols=4,
                as_html=False
            )
        except Exception as e:
            pytest.fail(f"display_alignment with SeqRecord list raised an exception: {e}")

    def test_seqrecord_with_max_seqs_limit(self):
        """Test that max_seqs works with SeqRecord lists"""
        mock_seqrecords = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("TGCA"), id="seq2"),
            SeqRecord(Seq("GGCC"), id="seq3"),
            SeqRecord(Seq("AATT"), id="seq4"),
        ]

        # Request only 2 sequences
        html_result = AlignmentViewer.get_alignment_html(
            mock_seqrecords,
            nseqs=2,
            ncols=4
        )

        # Should only contain the first 2 sequences
        assert "seq1" in html_result
        assert "seq2" in html_result
        assert "seq3" not in html_result
        assert "seq4" not in html_result

    def test_seqrecord_with_config_object(self):
        """Test using DisplayConfig object with SeqRecord lists"""
        mock_seqrecords = [
            SeqRecord(Seq("ACGTACGT"), id="sequence_1"),
            SeqRecord(Seq("TGCATGCA"), id="sequence_2"),
        ]

        config = DisplayConfig(
            nseqs=2,
            ncols=8,
            as_html=True,
            show_consensus=True
        )

        html_result = AlignmentViewer.get_alignment_html(mock_seqrecords, config=config)

        assert isinstance(html_result, str)
        assert "sequence_1" in html_result
        assert "sequence_2" in html_result
        # Check that individual nucleotides are present (they're in separate spans)
        assert ">A</span>" in html_result
        assert ">T</span>" in html_result

    def test_empty_seqrecord_list(self):
        """Test with empty list of SeqRecord objects"""
        empty_seqrecords = []

        with pytest.raises(Exception):  # Should raise an error for empty alignment
            AlignmentViewer.get_alignment_html(empty_seqrecords)

    @patch('Bio.SeqIO.parse')
    def test_real_seqio_parse_simulation(self, mock_seqio_parse):
        """Simulate the exact workflow: sequences = list(SeqIO.parse(...))"""
        # Mock what SeqIO.parse would return
        mock_seqio_parse.return_value = [
            SeqRecord(Seq("ACGTACGTACGT"), id="seq1", description="First sequence"),
            SeqRecord(Seq("TGCATGCATGCA"), id="seq2", description="Second sequence"),
            SeqRecord(Seq("GGCCGGCCGGCC"), id="seq3", description="Third sequence"),
        ]

        # Simulate the user's workflow
        with patch("builtins.open", mock_open(read_data="")):
            sequences = list(SeqIO.parse("../data/alignment.fasta", "fasta"))

            # This is the exact call from the notebook
            html_result = AlignmentViewer.get_alignment_html(
                sequences,
                nseqs=5,
                ncols=60,
                as_html=True
            )

            assert isinstance(html_result, str)
            assert "seq1" in html_result
            assert "seq2" in html_result
            assert "seq3" in html_result

    def test_seqrecord_preserves_original_functionality(self):
        """Ensure that adding SeqRecord support doesn't break file path input"""
        mock_fasta_content = ">seq1\nACGT\n>seq2\nTGCA\n"

        with patch("builtins.open", mock_open(read_data=mock_fasta_content)):
            with patch("Bio.SeqIO.parse") as mock_parse:
                mock_parse.return_value = [
                    SeqRecord(Seq("ACGT"), id="seq1"),
                    SeqRecord(Seq("TGCA"), id="seq2")
                ]

                # Test that file path input still works
                html_result = AlignmentViewer.get_alignment_html(
                    "test.fasta",
                    nseqs=2,
                    ncols=4
                )

                assert isinstance(html_result, str)
                assert "seq1" in html_result
                assert "seq2" in html_result
