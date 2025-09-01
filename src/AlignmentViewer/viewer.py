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
    def _prepare_alignment(alignment: Union[str, Path, TextIO, List],
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
        sequences = SequenceReader.parse(alignment, base_config.nseqs)

        # set limit for ncols:
        if base_config.ncols > len(sequences[0].sequence) or base_config.ncols == 0:
            base_config.ncols = len(sequences[0].sequence)

        # Calculate layout
        max_header_len = max(len(seq.header) for seq in sequences)

        return sequences, base_config, max_header_len

    @staticmethod
    def display_alignment(alignment: Union[str, Path, TextIO, List],
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
    def get_alignment_html(alignment: Union[str, Path, TextIO, List],
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
        consensus_css = ""
        if config.show_consensus:
            padding = ' ' * (max_header_len + 1)
            consensus_css, _ = self._generate_consensus_html(
                sequences, padding, config.ncols, config.start_pos,
                config.block_size, config.consensus_height, config.consensus_ignore_gaps
            )

        css += consensus_css + """
        </style>
        """

        html_parts.append(css)
        html_parts.append('<div class="alignment-container">')

        # Add ruler if requested
        padding = ' ' * (max_header_len + 1)

        if config.show_ruler:
            ruler = self._create_ruler(padding, config.ncols, config.start_pos, config.block_size)
            html_parts.append(f"<div class='sequence-row'>{ruler}</div>")

        # Add consensus bar chart if requested
        if config.show_consensus:
            _, consensus_html = self._generate_consensus_html(
                sequences, padding, config.ncols, config.start_pos,
                config.block_size, config.consensus_height, config.consensus_ignore_gaps
            )
            html_parts.append(consensus_html)

        # Add sequences
        consensus_seq = None
        if config.color_snps_only:
            calculator = ConsensusCalculator(sequences)
            consensus_seq = calculator.get_consensus_sequence(config.consensus_ignore_gaps)
        for sequence in sequences:
            colored_seq = self.formatter.format_sequence_html(
                sequence.sequence,
                config.ncols,
                config.block_size,
                config.start_pos,
                consensus=consensus_seq,
                color_snps_only=config.color_snps_only
            )
            html_parts.append(
                f"<div class='sequence-row'>"
                f"{sequence.header:<{max_header_len}} {colored_seq}"
                f"</div>"
            )

        html_parts.append('</div>')

        return ''.join(html_parts)

    def _generate_consensus_html(self, sequences: List[Sequence], padding: str,
                            ncols: int, start_pos: int, block_size: int,
                            consensus_height: str, consensus_ignore_gaps: bool) -> tuple[str, str]:
        """
        Generate HTML and CSS for consensus bar chart

        Returns:
            tuple: (css_string, html_string)
        """
        calculator = ConsensusCalculator(sequences)

        # Get consensus data for the visible range
        consensus_data = calculator.generate_consensus_bar_data(
            start_pos=start_pos,
            ncols=ncols,
            block_size=block_size,
            ignore_gaps=consensus_ignore_gaps
        )

        # Generate CSS for consensus with custom tooltip
        css = f"""
            .consensus-container {{
                height: {consensus_height};
                margin: 5px 0;
                border-bottom: 1px solid #ccc;
                display: flex;
                align-items: end;
                padding-bottom: 2px;
            }}
            .consensus-bar {{
                background-color: grey;
                margin: 0;
                padding: 0;
                box-sizing: border-box;
                position: relative;
                cursor: pointer;
                width: 1ch;
                font-family monospace;
                color: transparent;
            }}
            .consensus-bar:hover {{
                background-color: darkgrey;
            }}
            .consensus-bar::after {{
                content: attr(data-tooltip);
                position: absolute;
                bottom: 100%;
                left: 50%;
                transform: translateX(-50%);
                background-color: #333;
                color: white;
                padding: 4px 8px;
                border-radius: 4px;
                font-size: 11px;
                white-space: nowrap;
                visibility: hidden;
                opacity: 0;
                transition: opacity 0.3s ease;
                z-index: 1000;
                margin-bottom: 5px;
            }}
            .consensus-bar::before {{
                content: '';
                position: absolute;
                bottom: 100%;
                left: 50%;
                transform: translateX(-50%);
                border: 4px solid transparent;
                border-top-color: #333;
                visibility: hidden;
                opacity: 0;
                transition: opacity 0.3s ease;
                z-index: 1000;
                margin-bottom: 1px;
            }}
            .consensus-bar:hover::after,
            .consensus-bar:hover::before {{
                visibility: visible;
                opacity: 1;
            }}
            .consensus-spacer {{
                width: 1ch;
                font-family: monospace;
                height: 100%;
                background-color: transparent;
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}"""

        # Generate bar chart HTML with proper alignment using non-breaking spaces
        padding_html = '&nbsp;' * len(padding)
        consensus_html = f'<div class="sequence-row"><div class="consensus-container">{padding_html}'

        # Create bars for each position with proper spacing to match sequence formatting
        for i, (pos, agreement) in enumerate(consensus_data):
            # Add space before blocks (except the first position in each block)
            if i > 0 and i % block_size == 0:
                consensus_html += '<div class="consensus-spacer">&nbsp;</div>'

            # Calculate bar height as percentage of container height
            height_percent = agreement  # agreement is already 0-100
            bar_html = (
                f'<div class="consensus-bar" '
                f'style="height: {height_percent}%;" '
                f'data-tooltip="Pos {pos}: {agreement:.1f}%">'
                f'&nbsp;'
                f'</div>'
            )
            consensus_html += bar_html

        consensus_html += '</div></div>'
        return css, consensus_html

    def _display_notebook(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in a notebook with scrollable container"""
        html_content = self._generate_html(sequences, config, max_header_len)
        display(HTML(html_content))

    def _display_terminal(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment in terminal"""
        padding = ' ' * (max_header_len + 1)

        if config.show_ruler:
            ruler = self._create_ruler(padding, config.ncols, config.start_pos, config.block_size, is_html=False)
            print(ruler)

        consensus_seq = None
        if config.color_snps_only:
            calculator = ConsensusCalculator(sequences)
            consensus_seq = calculator.get_consensus_sequence(config.consensus_ignore_gaps)
        for sequence in sequences:
            colored_seq = self.formatter.format_sequence(
                sequence.sequence,
                config.ncols,
                config.block_size,
                config.start_pos,
                consensus=consensus_seq,
                color_snps_only=config.color_snps_only
            )
            print(f"{sequence.header:<{max_header_len}} {colored_seq}")

    def _create_ruler(self, padding: str, ncols: int, start_pos: int, block_size: int, is_html: bool=True) -> str:
        """Create a ruler string with column numbers and spaces"""
        ticks = ''
        numbers = []
        pos = start_pos

        # Create number line and tick marks simultaneously
        for i in range(ncols):
            if i % block_size == 0:
                num = str(pos + i)
                space = ' ' * (block_size - len(num))
                space += ' ' if i > 0 else ''
                numbers.append(num + space)
                ticks += '|'
                if i > 0:
                    ticks += ' '
            else:
                ticks += '-'
        whitespace = '&nbsp;' if is_html else ' '
        padding = whitespace * len(padding)
        return padding + ''.join(numbers) + '\n' + padding + ticks
