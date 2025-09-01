import pytest
from AlignmentViewer.sequence_formatter import SequenceFormatter
from AlignmentViewer.color_scheme import ColorScheme
from AlignmentViewer.consensus import ConsensusCalculator
from AlignmentViewer.sequence import Sequence

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


def test_snp_coloring_consensus_with_N():
    # Test case where consensus is all 'N' - no bases should be highlighted
    # Because N means "unknown", we shouldn't highlight any nucleotides as SNPs
    sequence = "ATCG"
    consensus = "NNNN"  # All N consensus means no valid consensus
    formatter = SequenceFormatter(ColorScheme.default())

    # Test HTML formatting
    result_html = formatter.format_sequence_html(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    # No bases should be colored since consensus is all N
    assert '<span' not in result_html
    assert result_html == "ATCG"

    # Test text formatting
    result_text = formatter.format_sequence(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    # No bases should be colored
    assert formatter.colors.nucleotides['A'] not in result_text
    assert formatter.colors.nucleotides['T'] not in result_text
    assert formatter.colors.nucleotides['C'] not in result_text
    assert formatter.colors.nucleotides['G'] not in result_text
    assert result_text == "ATCG"


def test_snp_coloring_mixed_N_consensus():
    # Test case where some positions have N consensus, others don't
    sequence = "ATCG"
    consensus = "ANNG"  # First and last positions have valid consensus
    formatter = SequenceFormatter(ColorScheme.default())

    # Test HTML formatting - only middle positions should NOT be colored (consensus is N)
    # Position 0: A vs A (no SNP), Position 1: T vs N (no SNP, N means unknown)
    # Position 2: C vs N (no SNP, N means unknown), Position 3: G vs G (no SNP)
    result_html = formatter.format_sequence_html(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    # No bases should be colored because positions 1,2 have N consensus (unknown) and positions 0,3 match
    assert '<span' not in result_html
    assert result_html == "ATCG"

    # Test text formatting
    result_text = formatter.format_sequence(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    # No bases should be colored
    assert formatter.colors.nucleotides['A'] not in result_text
    assert formatter.colors.nucleotides['T'] not in result_text
    assert formatter.colors.nucleotides['C'] not in result_text
    assert formatter.colors.nucleotides['G'] not in result_text
    assert result_text == "ATCG"


def test_snp_coloring_mixed_N_consensus_with_snps():
    # Test case where N consensus doesn't highlight, but real mismatches do
    sequence = "TTCG"
    consensus = "ANNG"  # Position 0: T vs A (SNP), Position 1: T vs N (no SNP),
                       # Position 2: C vs N (no SNP), Position 3: G vs G (no SNP)
    formatter = SequenceFormatter(ColorScheme.default())

    # Test HTML formatting - only first position should be colored (T vs A)
    result_html = formatter.format_sequence_html(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    span_count = result_html.count('<span')
    assert span_count == 1  # Only first T should be colored
    assert 'T</span>' in result_html

    # Test text formatting
    result_text = formatter.format_sequence(sequence, ncols=4, block_size=4, consensus=consensus, color_snps_only=True)
    # Only first T should be colored
    colored_T = formatter.colors.nucleotides['T'] + 'T' + formatter.colors.reset
    assert colored_T in result_text
    # Should have exactly one colored T
    assert result_text.count(formatter.colors.nucleotides['T']) == 1


def test_consensus_calculation_prioritizes_nucleotides_over_N():
    # Test that consensus calculation chooses real nucleotides over N
    sequences = [
        Sequence("seq1", "ANNN"),  # Position 0: A, others N
        Sequence("seq2", "ANNN"),  # Position 0: A, others N
        Sequence("seq3", "TNNN"),  # Position 0: T, others N
        Sequence("seq4", "ANNN"),  # Position 0: A, others N
    ]

    calc = ConsensusCalculator(sequences)

    # Position 0: A appears 3 times, T appears 1 time -> consensus should be A
    pos0_consensus = calc.calculate_position_consensus(0)
    assert pos0_consensus['most_common'] == 'A'

    # Position 1: All N -> consensus should be N (since no valid nucleotides)
    pos1_consensus = calc.calculate_position_consensus(1)
    assert pos1_consensus['most_common'] == 'N'

    # Test full consensus sequence
    consensus_seq = calc.get_consensus_sequence()
    assert consensus_seq == "ANNN"


def test_consensus_calculation_mixed_nucleotides_and_N():
    # Test that consensus ignores N when calculating most common nucleotide
    sequences = [
        Sequence("seq1", "ATCG"),
        Sequence("seq2", "NNCG"),  # N's at positions 0,1
        Sequence("seq3", "AACG"),
        Sequence("seq4", "ANCG"),  # N at position 1
    ]

    calc = ConsensusCalculator(sequences)

    # Position 0: A appears 3 times, ignore N -> consensus should be A
    pos0_consensus = calc.calculate_position_consensus(0)
    assert pos0_consensus['most_common'] == 'A'

    # Position 1: T appears 1 time, A appears 1 time, ignore 2 N's -> should pick most common among T,A
    # Since T and A each appear once, it should pick one of them (likely T as it comes first)
    pos1_consensus = calc.calculate_position_consensus(1)
    assert pos1_consensus['most_common'] in ['T', 'A']  # Either is valid since tied

    # Position 2: C appears 4 times -> consensus should be C
    pos2_consensus = calc.calculate_position_consensus(2)
    assert pos2_consensus['most_common'] == 'C'

