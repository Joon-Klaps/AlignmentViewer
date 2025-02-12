from typing import List, Union, TextIO, Optional
import sys
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

    def _display_notebook(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in a notebook with scrollable container"""
        html_parts = []

        # CSS for the container
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
        """Display alignment in terminal with width detection"""
        # Adjust ncols if it exceeds terminal width
        available_width = self.terminal_width - max_header_len - 2  # -2 for padding
        if config.ncols == 0 or config.ncols > available_width:
            config.ncols = available_width

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

    def _create_ruler(self, padding: str, width: int, start_pos: int = 0, step: int = 10) -> str:
        """Create a ruler string with column numbers and spaces using improved tick marks"""
        def calc_tick_indices(start, end, distance, min_distance):
            """Calculate positions for tick marks"""
            first_even_pos = (start // distance + 1) * distance
            if first_even_pos - start < min_distance:
                first_even_pos += distance
            return range(first_even_pos, end, distance)

        def make_one_tick(position, space, tickmark):
            """Create a single tick mark with position number"""
            return '{0:>{width}}{tickmark}'.format(position, width=space-1, tickmark=tickmark)

        # Use UTF-8 arrow if supported, else fallback to |
        tickmark = '↑' if sys.stdout.encoding == 'UTF-8' else '|'

        # Initial space for left margin
        index_bar = ' ' * (len(padding) - step + 1)

        # Add first column index
        index_bar += make_one_tick(start_pos, min(len(padding)+1, step), tickmark)

        # Add remaining ticks
        last_pos = start_pos
        for pos in calc_tick_indices(start_pos, start_pos + width, step, step):
            spacer = pos - last_pos - 1
            index_block = '{0:>{width}}{tickmark}'.format(pos, width=spacer, tickmark=tickmark)
            index_bar += index_block
            last_pos = pos

        return index_bar
