from typing import List, Union, TextIO, Optional
from pathlib import Path
from IPython.display import HTML, display

from .color_scheme import ColorScheme
from .config import DisplayConfig
from .sequence import Sequence
from .sequence_reader import SequenceReader
from .sequence_formatter import SequenceFormatter

class AlignmentViewer:
    def __init__(self, color_scheme: Optional[ColorScheme] = None):
        self.color_scheme = color_scheme or ColorScheme.default()
        self.formatter = SequenceFormatter(self.color_scheme)

    @staticmethod
    def display_alignment(alignment: Union[str, Path, TextIO, List[Sequence]],
                         config: Optional[DisplayConfig] = None,
                         **kwargs) -> None:
        """Display a colored alignment with optional configuration"""
        viewer = AlignmentViewer()

        # Create base config
        base_config = config or DisplayConfig()

        # Update config with any provided kwargs
        for key, value in kwargs.items():
            if hasattr(base_config, key):
                setattr(base_config, key, value)
            else:
                raise ValueError(f"Unknown configuration parameter: {key} \nPossible values are: {', '.join(base_config.__dict__.keys())}")

        # Validate config
        base_config.validate()

        # Parse sequences
        sequences = (alignment if isinstance(alignment, list)
                    else SequenceReader.parse(alignment, base_config.nseqs))

        # set limit for ncols:
        if base_config.ncols > len(sequences[0].sequence):
            base_config.ncols = 0 # set to all columns

        # Calculate layout
        max_header_len = max(len(seq.header) for seq in sequences)

        # Use config's as_html if set, otherwise detect
        as_html = (base_config.as_html if base_config.as_html is not None
                      else viewer._check_notebook())

        if as_html:
            viewer._display_notebook(sequences, base_config, max_header_len)
        else:
            viewer._display_terminal(sequences, base_config, max_header_len)

    def _check_notebook(self) -> bool:
        """Check if we're running in a Jupyter notebook"""
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True
            return False
        except NameError:
            return False

    def _calculate_display_width(self, sequence_length: int, block_size: int) -> int:
        """Calculate the actual display width considering spaces between blocks"""
        if sequence_length == 0:
            return 0
        # Add 1 space for every complete block (except the last one)
        num_blocks = (sequence_length - 1) // block_size
        return sequence_length + num_blocks

    def _create_ruler(self, padding: str, width: int, start_pos: int = 0, step: int = 10) -> str:
        """Create a ruler string with column numbers and spaces"""
        if width == 0:
            return ""

        ticks = []
        numbers = []
        pos = start_pos
        display_width = self._calculate_display_width(width, step)

        # Create the number line
        current_pos = 0
        while current_pos < display_width:
            if current_pos % (step + 1) == 0:  # +1 accounts for the space after each block
                num = str(pos + (current_pos * step // (step + 1)))
                numbers.append(num.ljust(step + 1))  # +1 for the space after block
            current_pos += 1

        # Create tick marks
        for i in range(display_width):
            if i % (step + 1) == 0:  # +1 accounts for the space after each block
                ticks.append('|')
            else:
                ticks.append('-')

        return ''.join(numbers).rstrip() + '\n' + padding + ''.join(ticks)

    def _display_notebook(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in a notebook with scrollable container"""
        html_parts = []

        # Add CSS with specific spacing adjustments
        css = f"""
        <style>
            .alignment-container {{
                font-family: monospace;
                white-space: pre;
                overflow-x: auto;
                max-height: {config.container_height};
                background-color: white;
                padding: 10px;
                border: 1px solid #ddd;
            }}
            .sequence-row {{
                line-height: 1.5;
                margin: 2px 0;
            }}
            .sequence-row span {{
                display: inline-block;
                width: 1ch;  /* Force consistent width */
            }}
        </style>
        """

        html_parts.append(css)
        html_parts.append('<div class="alignment-container">')

        # Add ruler if requested
        padding = ' ' * (max_header_len + 1)
        if config.show_ruler:
            ruler = self._create_ruler(padding, config.ncols, config.start_pos, config.block_size)
            html_parts.append(f"<div class='sequence-row'>{padding}{ruler}</div>")

        # Add sequences
        for sequence in sequences:
            colored_seq = self.formatter.format_sequence_html(
                sequence.sequence,
                config.ncols,
                config.block_size,
                config.start_pos
            )
            html_parts.append(
                f"<div class='sequence-row'>"
                f"{sequence.header:<{max_header_len}} {colored_seq}"
                f"</div>"
            )

        html_parts.append('</div>')

        # Display the HTML
        display(HTML(''.join(html_parts)))

    def _display_terminal(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in terminal"""
        padding = ' ' * (max_header_len + 1)

        if config.show_ruler:
            print(f"{padding}{self._create_ruler(padding, config.ncols, config.start_pos, config.block_size)}")

        for sequence in sequences:
            colored_seq = self.formatter.format_sequence(
                sequence.sequence,
                config.ncols,
                config.block_size,
                config.start_pos
            )
            print(f"{sequence.header:<{max_header_len}} {colored_seq}")
