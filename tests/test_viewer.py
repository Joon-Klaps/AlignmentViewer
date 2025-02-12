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