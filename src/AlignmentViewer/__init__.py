from .viewer import AlignmentViewer, ColorScheme, DisplayConfig, SequenceReader

__all__ = ['AlignmentViewer', 'ColorScheme', 'DisplayConfig', 'SequenceReader']

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"
