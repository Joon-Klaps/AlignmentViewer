from .viewer import AlignmentViewer, ColorScheme, DisplayConfig

__all__ = ['AlignmentViewer', 'ColorScheme', 'DisplayConfig']

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"
