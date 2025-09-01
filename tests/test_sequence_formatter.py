import pytest
from AlignmentViewer.sequence_formatter import SequenceFormatter
from AlignmentViewer.color_scheme import ColorScheme


def test_format_sequence():
    formatter = SequenceFormatter(ColorScheme.default())
    sequence = "ACGT"
    result = formatter.format_sequence(sequence, ncols=4, block_size=2)
    # Should contain ANSI color codes
    assert '\033[' in result
    assert 'A' in result
    assert 'C' in result


def test_format_sequence_html():
    formatter = SequenceFormatter(ColorScheme.default())
    sequence = "ACGT"
    result = formatter.format_sequence_html(sequence, ncols=4, block_size=2)
    # Should contain HTML spans
    assert '<span' in result
    assert 'style=' in result
    assert 'A' in result
    assert 'C' in result
