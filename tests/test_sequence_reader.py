import pytest
from unittest.mock import patch, mock_open
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from AlignmentViewer.sequence_reader import SequenceReader
from AlignmentViewer.sequence import Sequence


class TestSequenceReader:
    """Test cases for the SequenceReader class"""

    def test_parse_file_path(self):
        """Test parsing from a file path"""
        mock_fasta_content = ">seq1\nACGT\n>seq2\nTGCA\n"

        with patch("builtins.open", mock_open(read_data=mock_fasta_content)):
            with patch("Bio.SeqIO.parse") as mock_parse:
                # Mock SeqIO.parse to return SeqRecord objects
                mock_parse.return_value = [
                    SeqRecord(Seq("ACGT"), id="seq1"),
                    SeqRecord(Seq("TGCA"), id="seq2")
                ]

                sequences = SequenceReader.parse("test.fasta", max_seqs=0)

                assert len(sequences) == 2
                assert all(isinstance(seq, Sequence) for seq in sequences)
                assert sequences[0].header == "seq1"
                assert sequences[0].sequence == "ACGT"
                assert sequences[1].header == "seq2"
                assert sequences[1].sequence == "TGCA"

    def test_parse_list_seqrecord_objects(self):
        """Test parsing a list of BioPython SeqRecord objects"""
        seqrecords = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("TGCA"), id="seq2"),
            SeqRecord(Seq("GGCC"), id="seq3")
        ]

        sequences = SequenceReader.parse(seqrecords, max_seqs=0)

        assert len(sequences) == 3
        assert all(isinstance(seq, Sequence) for seq in sequences)
        assert sequences[0].header == "seq1"
        assert sequences[0].sequence == "ACGT"
        assert sequences[1].header == "seq2"
        assert sequences[1].sequence == "TGCA"
        assert sequences[2].header == "seq3"
        assert sequences[2].sequence == "GGCC"

    def test_parse_list_seqrecord_objects_with_limit(self):
        """Test parsing a list of BioPython SeqRecord objects with max_seqs limit"""
        seqrecords = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("TGCA"), id="seq2"),
            SeqRecord(Seq("GGCC"), id="seq3")
        ]

        sequences = SequenceReader.parse(seqrecords, max_seqs=2)

        assert len(sequences) == 2
        assert sequences[0].header == "seq1"
        assert sequences[1].header == "seq2"

    def test_parse_list_sequence_objects(self):
        """Test parsing a list of Sequence objects (should pass through)"""
        sequence_objects = [
            Sequence("seq1", "ACGT"),
            Sequence("seq2", "TGCA")
        ]

        sequences = SequenceReader.parse(sequence_objects, max_seqs=0)

        assert len(sequences) == 2
        assert sequences is sequence_objects  # Should be the same objects
        assert sequences[0].header == "seq1"
        assert sequences[0].sequence == "ACGT"

    def test_parse_list_sequence_objects_with_limit(self):
        """Test parsing a list of Sequence objects with max_seqs limit"""
        sequence_objects = [
            Sequence("seq1", "ACGT"),
            Sequence("seq2", "TGCA"),
            Sequence("seq3", "GGCC")
        ]

        sequences = SequenceReader.parse(sequence_objects, max_seqs=2)

        assert len(sequences) == 2
        assert sequences[0].header == "seq1"
        assert sequences[1].header == "seq2"

    def test_parse_empty_list(self):
        """Test parsing an empty list"""
        sequences = SequenceReader.parse([], max_seqs=0)
        assert len(sequences) == 0

    def test_parse_list_unsupported_objects(self):
        """Test parsing a list with unsupported object types"""
        unsupported_objects = ["string1", "string2"]

        with pytest.raises(ValueError, match="Unsupported list item type"):
            SequenceReader.parse(unsupported_objects, max_seqs=0)

    def test_parse_list_mixed_objects(self):
        """Test that the function properly detects object type from first item"""
        # Create a mix but it should only look at the first item
        mixed_objects = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            "not_a_sequence"  # This would cause issues but should be handled by from_seqrecord
        ]

        with pytest.raises(Exception):  # from_seqrecord will fail on the string
            SequenceReader.parse(mixed_objects, max_seqs=0)

    def test_guess_format(self):
        """Test format guessing from file extensions"""
        assert SequenceReader._guess_format("test.fasta") == "fasta"
        assert SequenceReader._guess_format("test.fa") == "fasta"
        assert SequenceReader._guess_format("test.sto") == "stockholm"
        assert SequenceReader._guess_format("test.gb") == "genbank"
        assert SequenceReader._guess_format("test.gbk") == "genbank"
        assert SequenceReader._guess_format("test.phylip") == "phylip"
        assert SequenceReader._guess_format("test.phy") == "phylip"
        assert SequenceReader._guess_format("test.clustal") == "clustal"
        assert SequenceReader._guess_format("test.aln") == "clustal"
        assert SequenceReader._guess_format("test.embl") == "embl"
        assert SequenceReader._guess_format("test.unknown") == "fasta"  # default

    def test_parse_file_with_max_seqs(self):
        """Test parsing file with max_seqs limit"""
        mock_fasta_content = ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nGGCC\n"

        with patch("builtins.open", mock_open(read_data=mock_fasta_content)):
            with patch("Bio.SeqIO.parse") as mock_parse:
                mock_parse.return_value = [
                    SeqRecord(Seq("ACGT"), id="seq1"),
                    SeqRecord(Seq("TGCA"), id="seq2"),
                    SeqRecord(Seq("GGCC"), id="seq3")
                ]

                sequences = SequenceReader.parse("test.fasta", max_seqs=2)

                assert len(sequences) == 2
                assert sequences[0].header == "seq1"
                assert sequences[1].header == "seq2"

    def test_parse_error_handling(self):
        """Test error handling for invalid files"""
        with patch("builtins.open", side_effect=FileNotFoundError("File not found")):
            with pytest.raises(ValueError, match="Error parsing sequence file"):
                SequenceReader.parse("nonexistent.fasta", max_seqs=0)

    def test_max_seqs_greater_than_available(self):
        """Test when max_seqs is greater than available sequences"""
        seqrecords = [
            SeqRecord(Seq("ACGT"), id="seq1"),
            SeqRecord(Seq("TGCA"), id="seq2")
        ]

        # Should return all sequences when max_seqs > available
        sequences = SequenceReader.parse(seqrecords, max_seqs=10)
        assert len(sequences) == 2

    def test_parse_seqrecord_with_description(self):
        """Test parsing SeqRecord with description"""
        seqrecord = SeqRecord(
            Seq("ACGT"),
            id="seq1",
            description="Test sequence"
        )

        sequences = SequenceReader.parse([seqrecord], max_seqs=0)
        assert len(sequences) == 1
        assert sequences[0].header == "seq1"  # Should use id, not description
        assert sequences[0].sequence == "ACGT"
