import pytest
from pathlib import Path

def test_import():
    """Test that we can import the package"""
    from AlignmentViewer import AlignmentViewer
    assert AlignmentViewer is not None

def test_import_consensus():
    """Test that we can import the ConsensusCalculator"""
    from AlignmentViewer import ConsensusCalculator
    assert ConsensusCalculator is not None

def test_display_alignment():
    """Test that we can display an alignment"""
    from AlignmentViewer import AlignmentViewer

    # Get the path to the test data
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"
    assert data_path.exists(), f"Test data not found at {data_path}"

    # Basic test that display doesn't raise an error
    try:
        AlignmentViewer.display_alignment(str(data_path), nseqs=3, ncols=50, as_html=False)
    except Exception as e:
        pytest.fail(f"display_alignment raised an exception: {e}")

def test_display_alignment_with_html():
    """Test that we can display an alignment with HTML"""
    from AlignmentViewer import AlignmentViewer

    # Get the path to the test data
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"
    assert data_path.exists(), f"Test data not found at {data_path}"

    # Basic test that display doesn't raise an error
    try:
        AlignmentViewer.display_alignment(str(data_path), as_html=True)
    except Exception as e:
        pytest.fail(f"display_alignment raised an exception: {e}")

def test_display_alignment_with_clustal_input():
    """Test that we can display an alignment with clustal input"""
    from AlignmentViewer import AlignmentViewer

    # Get the path to the test data
    data_path = Path(__file__).parent.parent / "data" / "alignment.aln"
    assert data_path.exists(), f"Test data not found at {data_path}"

    # Basic test that display doesn't raise an error
    try:
        AlignmentViewer.display_alignment(str(data_path), as_html=False)
    except Exception as e:
        pytest.fail(f"display_alignment raised an exception: {e}")

def test_fail_nseqs_negative():
    """Test that display_alignment fails with negative nseqs"""
    from AlignmentViewer import AlignmentViewer
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"
    with pytest.raises(ValueError):
        AlignmentViewer.display_alignment(str(data_path), nseqs=-1, as_html=False)

def test_fail_nseqs_float():
    """Test that display_alignment fails with float nseqs"""
    from AlignmentViewer import AlignmentViewer
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"

    with pytest.raises(ValueError):
        AlignmentViewer.display_alignment(str(data_path), nseqs=1.5, as_html=False)

def test_fail_ncols_negative():
    """Test that display_alignment fails with negative ncols"""
    from AlignmentViewer import AlignmentViewer
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"
    with pytest.raises(ValueError):
        AlignmentViewer.display_alignment(str(data_path), ncols=-1, as_html=False)

def test_fail_ncols_too_large():
    """Test that display_alignment fails with extremely large ncols"""
    from AlignmentViewer import AlignmentViewer
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"

    AlignmentViewer.display_alignment(str(data_path), ncols=50000000000, as_html=False)

def test_fail_startpos_negative():
    """Test that display_alignment fails with negative start_pos"""
    from AlignmentViewer import AlignmentViewer
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"
    with pytest.raises(ValueError):
        AlignmentViewer.display_alignment(str(data_path), start_pos=-1, as_html=False)

def test_fail_startpos_float():
    """Test that display_alignment fails with float start_pos"""
    from AlignmentViewer import AlignmentViewer
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"
    with pytest.raises(TypeError):
        AlignmentViewer.display_alignment(str(data_path), start_pos=0.5, as_html=False)

def test_get_alignment_html():
    """Test that we can get HTML output from alignment"""
    from AlignmentViewer import AlignmentViewer
    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"

    html_output = AlignmentViewer.get_alignment_html(str(data_path), nseqs=3, ncols=50)

    assert isinstance(html_output, str)
    assert len(html_output) > 0
    assert "alignment-container" in html_output
    assert "sequence-row" in html_output

def test_consensus_calculator_basic():
    """Test basic ConsensusCalculator functionality"""
    from AlignmentViewer import ConsensusCalculator, Sequence

    # Create test sequences with known consensus
    sequences = [
        Sequence("seq1", "ATCG"),
        Sequence("seq2", "ATCG"),
        Sequence("seq3", "ATCG"),
        Sequence("seq4", "ACCG"),  # Different at position 1
    ]

    calculator = ConsensusCalculator(sequences)

    # Test position 0 (should be 100% agreement - all A)
    pos_0_consensus = calculator.calculate_position_consensus(0)
    assert pos_0_consensus['most_common'] == 'A'
    assert pos_0_consensus['agreement_percentage'] == 100.0

    # Test position 1 (should be 75% agreement - 3 T, 1 C)
    pos_1_consensus = calculator.calculate_position_consensus(1)
    assert pos_1_consensus['most_common'] == 'T'
    assert pos_1_consensus['agreement_percentage'] == 75.0

    # Test agreement percentages for all positions
    agreements = calculator.get_consensus_agreement_percentages()
    assert len(agreements) == 4
    assert agreements[0] == 100.0  # All A
    assert agreements[1] == 75.0   # 3 T, 1 C
    assert agreements[2] == 100.0  # All C
    assert agreements[3] == 100.0  # All G

def test_consensus_calculator_with_gaps():
    """Test ConsensusCalculator with gap characters"""
    from AlignmentViewer import ConsensusCalculator, Sequence

    sequences = [
        Sequence("seq1", "ATG-"),
        Sequence("seq2", "ATG-"),
        Sequence("seq3", "AT--"),
        Sequence("seq4", "AT--"),
    ]

    calculator = ConsensusCalculator(sequences)

    # Test with gaps ignored (default)
    pos_2_consensus = calculator.calculate_position_consensus(2, ignore_gaps=True)
    assert pos_2_consensus['most_common'] == 'G'
    assert pos_2_consensus['agreement_percentage'] == 100.0
    assert pos_2_consensus['non_gap_sequences'] == 2

    # Test with gaps included
    pos_2_consensus_with_gaps = calculator.calculate_position_consensus(2, ignore_gaps=False)
    assert pos_2_consensus_with_gaps['most_common'] in ['G', '-']  # Could be either depending on tie-breaking
    assert pos_2_consensus_with_gaps['non_gap_sequences'] == 4

def test_consensus_calculator_validation():
    """Test ConsensusCalculator input validation"""
    from AlignmentViewer import ConsensusCalculator, Sequence

    # Test empty sequences
    with pytest.raises(ValueError, match="No sequences provided"):
        ConsensusCalculator([])

    # Test sequences with different lengths
    sequences = [
        Sequence("seq1", "ATCG"),
        Sequence("seq2", "ATCGAA"),  # Different length
    ]
    with pytest.raises(ValueError, match="All sequences must have the same length"):
        ConsensusCalculator(sequences)

    # Test position out of range
    valid_sequences = [
        Sequence("seq1", "ATCG"),
        Sequence("seq2", "ATCG"),
    ]
    calculator = ConsensusCalculator(valid_sequences)

    with pytest.raises(ValueError, match="Position .* out of range"):
        calculator.calculate_position_consensus(-1)

    with pytest.raises(ValueError, match="Position .* out of range"):
        calculator.calculate_position_consensus(10)

def test_consensus_sequence_generation():
    """Test consensus sequence generation"""
    from AlignmentViewer import ConsensusCalculator, Sequence

    sequences = [
        Sequence("seq1", "ATCG"),
        Sequence("seq2", "ATCG"),
        Sequence("seq3", "ACCG"),  # Different at position 1
    ]

    calculator = ConsensusCalculator(sequences)
    consensus_seq = calculator.get_consensus_sequence()

    assert consensus_seq == "ATCG"  # T wins at position 1 (2 vs 1)

def test_consensus_bar_data():
    """Test consensus bar data generation for visualization"""
    from AlignmentViewer import ConsensusCalculator, Sequence

    sequences = [
        Sequence("seq1", "ATCGATCG"),
        Sequence("seq2", "ATCGATCG"),
        Sequence("seq3", "ACCGATCG"),  # Different at position 1
    ]

    calculator = ConsensusCalculator(sequences)

    # Test full range
    bar_data = calculator.generate_consensus_bar_data()
    assert len(bar_data) == 8

    # Check specific values
    assert bar_data[0] == (0, 100.0)  # Position 0, 100% agreement
    assert bar_data[1][0] == 1  # Position 1
    assert abs(bar_data[1][1] - 66.67) < 0.1  # ~66.67% agreement (2/3)

    # Test partial range
    bar_data_partial = calculator.generate_consensus_bar_data(start_pos=2, ncols=3)
    assert len(bar_data_partial) == 3
    assert bar_data_partial[0][0] == 2  # Starting at position 2

def test_display_alignment_with_consensus():
    """Test that we can display an alignment with consensus"""
    from AlignmentViewer import AlignmentViewer

    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"
    assert data_path.exists(), f"Test data not found at {data_path}"

    # Test that display with consensus doesn't raise an error
    try:
        AlignmentViewer.display_alignment(
            str(data_path),
            nseqs=3,
            ncols=50,
            show_consensus=True,
            consensus_height="60px",
            as_html=False
        )
    except Exception as e:
        pytest.fail(f"display_alignment with consensus raised an exception: {e}")

def test_get_alignment_html_with_consensus():
    """Test HTML generation with consensus visualization"""
    from AlignmentViewer import AlignmentViewer

    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"

    # Test HTML generation with consensus
    html_output = AlignmentViewer.get_alignment_html(
        str(data_path),
        nseqs=3,
        ncols=50,
        show_consensus=True,
        consensus_height="60px",
        consensus_ignore_gaps=True
    )

    assert isinstance(html_output, str)
    assert len(html_output) > 0

    # Check that consensus-specific elements are present
    assert "consensus-container" in html_output
    assert "consensus-bar" in html_output
    assert "60px" in html_output  # Custom height

    # Check CSS classes
    assert ".consensus-container" in html_output
    assert ".consensus-bar" in html_output

def test_consensus_configuration_options():
    """Test various consensus configuration options"""
    from AlignmentViewer import AlignmentViewer, DisplayConfig

    data_path = Path(__file__).parent.parent / "data" / "alignment.fasta"

    # Test different consensus configurations
    configs = [
        DisplayConfig(show_consensus=True, consensus_height="40px", consensus_ignore_gaps=True),
        DisplayConfig(show_consensus=True, consensus_height="80px", consensus_ignore_gaps=False),
        DisplayConfig(show_consensus=False),  # Disabled consensus
    ]

    for config in configs:
        html_output = AlignmentViewer.get_alignment_html(str(data_path), config, nseqs=3)

        if config.show_consensus:
            assert "consensus-container" in html_output
            assert config.consensus_height in html_output
        else:
            assert "consensus-container" not in html_output

def test_consensus_with_sequence_objects():
    """Test consensus functionality with Sequence objects directly"""
    from AlignmentViewer import AlignmentViewer, Sequence

    # Create test sequences
    sequences = [
        Sequence("seq1", "ATCGATCG"),
        Sequence("seq2", "ATCGATCG"),
        Sequence("seq3", "ACCGATCG"),  # Different at position 1
        Sequence("seq4", "ACCGATCG"),
    ]

    # Test HTML generation
    html_output = AlignmentViewer.get_alignment_html(
        sequences,
        show_consensus=True,
        consensus_height="50px"
    )

    assert "consensus-container" in html_output
    assert "consensus-bar" in html_output

    # Test display
    try:
        AlignmentViewer.display_alignment(
            sequences,
            show_consensus=True,
            as_html=False
        )
    except Exception as e:
        pytest.fail(f"display_alignment with sequence objects and consensus raised an exception: {e}")
