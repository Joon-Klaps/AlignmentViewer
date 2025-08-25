import pytest
import numpy as np
from AlignmentViewer import AlignmentViewer
from AlignmentViewer.config import DisplayConfig
from AlignmentViewer.sequence import Sequence


class TestPlotlyIntegration:
    """Test suite for Plotly integration features"""

    def test_config_as_plotly_option(self):
        """Test that DisplayConfig accepts as_plotly parameter"""
        config = DisplayConfig(as_plotly=True)
        assert config.as_plotly is True

    def test_config_warning_both_html_and_plotly(self):
        """Test warning when both as_html and as_plotly are True"""
        config = DisplayConfig(as_html=True, as_plotly=True)
        with pytest.warns(UserWarning, match="Both as_html and as_plotly are set to True"):
            config.validate()


    def test_create_plotly_figure_method_exists(self):
        """Test that _create_plotly_figure method exists in AlignmentViewer"""
        viewer = AlignmentViewer()
        assert hasattr(viewer, '_create_plotly_figure')
        assert callable(getattr(viewer, '_create_plotly_figure'))

    def test_plotly_figure_creation_with_new_implementation(self, monkeypatch):
        """Test that figure creation works with new discrete color implementation"""
        # Mock all plotly components for new implementation
        class DummyHeatmap:
            def __init__(self, *args, **kwargs):
                self.z = kwargs.get('z', [])
                self.x = kwargs.get('x', [])
                self.y = kwargs.get('y', [])
                self.text = kwargs.get('text', [])

        class DummyFig:
            def __init__(self, data=None):
                self.data = data or []
            def show(self):
                return True
            def update_layout(self, *args, **kwargs):
                return self
            def add_trace(self, *args, **kwargs):
                return self

        class DummyBar:
            def __init__(self, *args, **kwargs):
                pass

        class DummySubplots:
            def add_trace(self, *args, **kwargs):
                return self
            def update_layout(self, *args, **kwargs):
                return self
            def update_xaxes(self, *args, **kwargs):
                return self
            def update_yaxes(self, *args, **kwargs):
                return self

        def mock_make_subplots(*args, **kwargs):
            return DummySubplots()

        monkeypatch.setattr("plotly.graph_objects.Figure", DummyFig)
        monkeypatch.setattr("plotly.graph_objects.Heatmap", DummyHeatmap)
        monkeypatch.setattr("plotly.graph_objects.Bar", DummyBar)
        monkeypatch.setattr("plotly.subplots.make_subplots", mock_make_subplots)

        sequences = [
            Sequence("seq1", "ACGT"),
            Sequence("seq2", "TGCA")
        ]

        viewer = AlignmentViewer()
        config = DisplayConfig(as_plotly=True)

        # Should create figure without errors using new implementation
        fig = viewer._create_plotly_figure(sequences, config)
        assert fig is not None

    def test_get_alignment_plotly_static_method(self, monkeypatch):
        """Test AlignmentViewer.get_alignment_plotly static method"""
        # Mock plotly imports
        class DummyFig:
            def show(self):
                return True
            def update_layout(self, *args, **kwargs):
                return self
            def to_html(self, include_plotlyjs=True):
                return "<html>Mock figure</html>"

        monkeypatch.setattr("plotly.express.imshow", lambda *args, **kwargs: DummyFig())

        sequences = [
            Sequence("seq1", "ACGT"),
            Sequence("seq2", "TGCA")
        ]

        fig = AlignmentViewer.get_alignment_plotly(sequences)
        assert fig is not None


    def test_display_alignment_with_as_plotly_config(self, monkeypatch):
        """Test display_alignment with as_plotly=True in config"""
        # Mock plotly imports and figure
        class DummyFig:
            def show(self):
                return True
            def update_layout(self, *args, **kwargs):
                return self

        monkeypatch.setattr("plotly.express.imshow", lambda *args, **kwargs: DummyFig())

        sequences = [
            Sequence("seq1", "ACGT"),
            Sequence("seq2", "TGCA")
        ]

        config = DisplayConfig(as_plotly=True)

        # Should complete without errors
        result = AlignmentViewer.display_alignment(sequences, config=config)
        assert result is None  # display_alignment returns None

    def test_consensus_bar_with_new_plotly_implementation(self, monkeypatch):
        """Test Plotly output with consensus bar enabled using new implementation"""
        # Mock plotly components for new implementation
        class DummyHeatmap:
            def __init__(self, *args, **kwargs):
                pass

        class DummyBar:
            def __init__(self, *args, **kwargs):
                pass

        class DummySubplots:
            def add_trace(self, *args, **kwargs):
                return self
            def update_layout(self, *args, **kwargs):
                return self
            def update_xaxes(self, *args, **kwargs):
                return self
            def update_yaxes(self, *args, **kwargs):
                return self
            def show(self):
                return True

        def mock_make_subplots(*args, **kwargs):
            return DummySubplots()

        monkeypatch.setattr("plotly.graph_objects.Heatmap", DummyHeatmap)
        monkeypatch.setattr("plotly.graph_objects.Bar", DummyBar)
        monkeypatch.setattr("plotly.subplots.make_subplots", mock_make_subplots)

        sequences = [
            Sequence("seq1", "ACGT"),
            Sequence("seq2", "TGCA"),
            Sequence("seq3", "ACGT"),
        ]

        config = DisplayConfig(as_plotly=True, show_consensus=True)

        # Should complete without errors
        result = AlignmentViewer.display_alignment(sequences, config=config)
        assert result is None
