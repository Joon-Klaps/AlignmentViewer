import pytest
from AlignmentViewer.sequence_formatter import SequenceFormatter
from AlignmentViewer.color_scheme import ColorScheme

def test_snp_coloring_html():
    # Consensus: ACGT
    sequence = "ACGA"
    consensus = "ACGT"
    formatter = SequenceFormatter(ColorScheme.default())
    # Only last base is a SNP
    result = formatter.format_sequence_html(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    # Should only color the last base (A), not the first three
    assert result.count('<span') == 1
    assert 'A' in result
    assert 'G' in result
    assert 'T' not in result  # T is consensus, not in sequence


def test_snp_coloring_text():
    # Consensus: AAAA
    sequence = "AATA"
    consensus = "AAAA"
    formatter = SequenceFormatter(ColorScheme.default())
    # Only third base is a SNP
    result = formatter.format_sequence(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    # Should only color the third base (T)
    colored_T = formatter.colors.nucleotides['T'] + 'T' + formatter.colors.reset
    assert colored_T in result
    # The other bases should not be colored
    assert result.count(formatter.colors.nucleotides['A']) == 0
    assert result.count(formatter.colors.nucleotides['T']) == 1
