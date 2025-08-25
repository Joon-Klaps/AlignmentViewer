from typing import List, Union, TextIO, Optional
from pathlib import Path
from IPython.display import HTML, display

from .color_scheme import ColorScheme
from .config import DisplayConfig
from .sequence import Sequence
from .sequence_reader import SequenceReader
from .sequence_formatter import SequenceFormatter
from .consensus import ConsensusCalculator

class AlignmentViewer:
    def __init__(self, color_scheme: Optional[ColorScheme] = None):
        self.color_scheme = color_scheme or ColorScheme.default()
        self.formatter = SequenceFormatter(self.color_scheme)

    @staticmethod
    def _prepare_alignment(alignment: Union[str, Path, TextIO, List[Sequence]],
                          config: Optional[DisplayConfig] = None,
                          **kwargs) -> tuple[List[Sequence], DisplayConfig, int]:
        """Shared preparation logic for alignment processing"""
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
        if base_config.ncols > len(sequences[0].sequence) or base_config.ncols == 0:
            base_config.ncols = len(sequences[0].sequence)

        # Calculate layout
        max_header_len = max(len(seq.header) for seq in sequences)

        return sequences, base_config, max_header_len

    @staticmethod
    def display_alignment(alignment: Union[str, Path, TextIO, List[Sequence]],
                         config: Optional[DisplayConfig] = None,
                         **kwargs) -> None:
        """Display a colored alignment with optional configuration"""
        viewer = AlignmentViewer()
        sequences, base_config, max_header_len = AlignmentViewer._prepare_alignment(
            alignment, config, **kwargs
        )

        # Use config's as_html if set, otherwise detect
        as_html = (base_config.as_html if base_config.as_html is not None
                      else viewer._check_notebook())

        if as_html:
            viewer._display_notebook(sequences, base_config, max_header_len)
        else:
            viewer._display_terminal(sequences, base_config, max_header_len)

    @staticmethod
    def get_alignment_html(alignment: Union[str, Path, TextIO, List[Sequence]],
                          config: Optional[DisplayConfig] = None,
                          **kwargs) -> str:
        """Generate raw HTML for a colored alignment that can be embedded in other platforms"""
        viewer = AlignmentViewer()
        sequences, base_config, max_header_len = AlignmentViewer._prepare_alignment(
            alignment, config, **kwargs
        )

        return viewer._generate_html(sequences, base_config, max_header_len)

    def _check_notebook(self) -> bool:
        """Check if we're running in a Jupyter notebook"""
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True
            return False
        except NameError:
            return False

    def _generate_html(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> str:
        """Generate raw HTML for alignment (shared by notebook display and external use)"""
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
            }}"""

        # Add consensus CSS only if consensus is enabled
        if config.show_consensus:
            css += f"""
            .consensus-container {{
                height: {config.consensus_height};
                margin: 5px 0;
                border-bottom: 1px solid #ccc;
                display: flex;
                align-items: end;
                padding-bottom: 2px;
            }}
            .consensus-bar {{
                width: 1ch;
                background-color: #4CAF50;
                margin-right: 0;
                border-right: 1px solid #fff;
            }}
            .consensus-bar:hover {{
                background-color: #45a049;
            }}"""

        css += """
        </style>
        """

        html_parts.append(css)
        html_parts.append('<div class="alignment-container">')

        # Add ruler if requested
        padding = ' ' * (max_header_len + 1)

        if config.show_ruler:
            ruler = self._create_ruler(padding, config.ncols, config.start_pos, config.block_size)
            html_parts.append(f"<div class='sequence-row'>{padding}{ruler}</div>")

        # Add consensus bar chart if requested
        if config.show_consensus:
            consensus_html = self._generate_consensus_html(sequences, config, max_header_len)
            html_parts.append(consensus_html)

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

        return ''.join(html_parts)

    def _generate_consensus_html(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> str:
        """Generate HTML for consensus bar chart"""
        calculator = ConsensusCalculator(sequences)

        # Get consensus data for the visible range
        consensus_data = calculator.generate_consensus_bar_data(
            start_pos=config.start_pos,
            ncols=config.ncols,
            block_size=config.block_size,
            ignore_gaps=config.consensus_ignore_gaps
        )

        # Generate bar chart HTML
        padding = ' ' * (max_header_len + 1)
        consensus_html = f'<div class="sequence-row">{padding}<div class="consensus-container">'

        # Create bars for each position
        for pos, agreement in consensus_data:
            # Calculate bar height as percentage of container height
            height_percent = agreement  # agreement is already 0-100
            bar_html = (
                f'<div class="consensus-bar" '
                f'style="height: {height_percent}%;" '
                f'title="Position {pos + 1}: {agreement:.1f}% agreement">'
                f'</div>'
            )
            consensus_html += bar_html

        consensus_html += '</div></div>'
        return consensus_html

    def _display_notebook(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in a notebook with scrollable container"""
        html_content = self._generate_html(sequences, config, max_header_len)
        display(HTML(html_content))

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

    def _create_ruler(self, padding: str, width: int, start_pos: int = 0, step: int = 10) -> str:
        """Create a ruler string with column numbers and spaces"""
        ticks = ''
        numbers = []
        pos = start_pos

        # Create number line and tick marks simultaneously
        for i in range(width):
            if i % step == 0:
                num = str(pos + i)
                space = ' ' * (step - len(num))
                space += ' ' if i > 0 else ''
                numbers.append(num + space)
                ticks += '|'
                if i > 0:
                    ticks += ' '
            else:
                ticks += '-'

        return ''.join(numbers) + '\n' + padding + ticks
