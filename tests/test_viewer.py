import pytest
from pathlib import Path

def test_import():
    """Test that we can import the package"""
    from AlignmentViewer import AlignmentViewer
    assert AlignmentViewer is not None

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
