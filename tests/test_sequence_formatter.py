import numpy as np
import pytest
from AlignmentViewer.sequence_formatter import SequenceFormatter
from AlignmentViewer.viewer import AlignmentViewer
from AlignmentViewer.color_scheme import ColorScheme
from AlignmentViewer.config import DisplayConfig
from AlignmentViewer.sequence import Sequence


def test_alignment_to_dataframe():
    sequences = [
        "ACGT",
        "TGCA",
        "CCCC",
        "GGGG"
    ]
    formatter = SequenceFormatter(ColorScheme.default())
    df = formatter.alignment_to_dataframe(sequences)

    # Test DataFrame properties
    assert df.shape == (4, 4)
    assert list(df.columns) == [1, 2, 3, 4]
    assert list(df.index) == ['Seq_1', 'Seq_2', 'Seq_3', 'Seq_4']
    assert df.iloc[0, 0] == 'A'
    assert df.iloc[1, 1] == 'G'


def test_alignment_to_dataframe_with_names():
    sequences = ["ACGT", "TGCA"]
    names = ["Sequence_A", "Sequence_B"]
    formatter = SequenceFormatter(ColorScheme.default())
    df = formatter.alignment_to_dataframe(sequences, names)

    assert list(df.index) == ["Sequence_A", "Sequence_B"]
    assert df.iloc[0, 0] == 'A'


def test_alignment_to_dataframe_raises_on_unequal_length():
    sequences = ["ACGT", "TGCA", "CC"]
    formatter = SequenceFormatter(ColorScheme.default())
    with pytest.raises(ValueError):
        formatter.alignment_to_dataframe(sequences)


def test_plotly_heatmap_integration(monkeypatch):
    # Mock plotly to avoid opening browser during test
    class DummyFig:
        def show(self):
            return True
        def update_layout(self, *args, **kwargs):
            return self

    monkeypatch.setattr("plotly.express.imshow", lambda *args, **kwargs: DummyFig())

    sequences = [
        Sequence("seq1", "ACGT"),
        Sequence("seq2", "TGCA"),
        Sequence("seq3", "CCCC"),
        Sequence("seq4", "GGGG")
    ]

    # Test integrated Plotly functionality
    result = AlignmentViewer.display_alignment(sequences, as_plotly=True)
    assert result is None  # display_alignment returns None


def test_config_as_plotly_validation():
    config = DisplayConfig(as_html=True, as_plotly=True)
    # Should warn but not raise error
    with pytest.warns(UserWarning):
        config.validate()


def test_get_alignment_plotly(monkeypatch):
    # Mock plotly imports
    class DummyFig:
        def show(self):
            return True
        def update_layout(self, *args, **kwargs):
            return self
        def to_html(self, include_plotlyjs=True):
            return "<html>Mock plotly figure</html>"

    def mock_imshow(*args, **kwargs):
        return DummyFig()

    monkeypatch.setattr("plotly.express.imshow", mock_imshow)

    # Create test sequences
    sequences = [
        Sequence("seq1", "ACGT"),
        Sequence("seq2", "TGCA")
    ]

    fig = AlignmentViewer.get_alignment_plotly(sequences)
    assert fig is not None

def test_display_alignment_with_plotly(monkeypatch):
    # Mock plotly imports and display
    class DummyFig:
        def show(self):
            return True
        def update_layout(self, *args, **kwargs):
            return self
        def to_html(self, include_plotlyjs=True):
            return "<html>Mock plotly figure</html>"

    def mock_imshow(*args, **kwargs):
        return DummyFig()

    monkeypatch.setattr("plotly.express.imshow", mock_imshow)

    # Create test sequences
    sequences = [
        Sequence("seq1", "ACGT"),
        Sequence("seq2", "TGCA")
    ]

    config = DisplayConfig(as_plotly=True)

    # Should not raise any errors
    AlignmentViewer.display_alignment(sequences, config=config)
