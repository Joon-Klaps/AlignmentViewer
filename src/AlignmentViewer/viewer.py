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

        # Integrated Plotly output
        as_plotly = getattr(base_config, 'as_plotly', False)
        if as_plotly:
            viewer._display_plotly(sequences, base_config, max_header_len)
            return

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

    @staticmethod
    def get_alignment_plotly(alignment: Union[str, Path, TextIO, List],
                            config: Optional[DisplayConfig] = None,
                            **kwargs):
        """Generate a Plotly figure object for a colored alignment"""
        viewer = AlignmentViewer()
        sequences, base_config, max_header_len = AlignmentViewer._prepare_alignment(
            alignment, config, **kwargs
        )

        return viewer._create_plotly_figure(sequences, base_config)

    @staticmethod
    def get_alignment_plotly_html(alignment: Union[str, Path, TextIO, List],
                                 config: Optional[DisplayConfig] = None,
                                 **kwargs) -> str:
        """Generate Plotly HTML string for a colored alignment"""
        fig = AlignmentViewer.get_alignment_plotly(alignment, config, **kwargs)
        return fig.to_html(include_plotlyjs=True)

    def _check_notebook(self) -> bool:
        """Check if we're running in a Jupyter notebook"""
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True
            return False
        except NameError:
            return False

    def _display_plotly(self, sequences: List[Sequence], config: DisplayConfig, max_header_len: int) -> None:
        """Display alignment using Plotly"""
        fig = self._create_plotly_figure(sequences, config)

        # Display locally (chart_studio requires authentication and can cause issues)
        fig.show()

    def _create_plotly_figure(self, sequences: List[Sequence], config: DisplayConfig):
        """Create Plotly figure for alignment visualization"""
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        # Convert sequences to DataFrame for better labeling
        sequence_names = [seq.header for seq in sequences]
        sequence_strings = [seq.sequence for seq in sequences]
        df = self.formatter.alignment_to_dataframe(sequence_strings, sequence_names)

        # Apply performance optimizations
        df_optimized = self._optimize_data_for_plotly(df, config)

        # Get color mapping from color scheme using the new method
        color_map = self.color_scheme.get_plotly_colors()

        # Create numerical matrix for plotting and color matrix
        positions = df_optimized.columns.tolist()
        seq_names = df_optimized.index.tolist()
        total_positions = len(positions)

        # Create arrays for heatmap
        z_values = []  # numerical values for heatmap
        colors = []    # color values
        text_values = []  # text to display

        for i, seq_name in enumerate(seq_names):
            row_z = []
            row_colors = []
            row_text = []
            for j, pos in enumerate(positions):
                nucleotide = df_optimized.iloc[i, j]
                # Map nucleotide to number for z-values
                nuc_num = self._nucleotide_to_number(nucleotide)
                row_z.append(nuc_num)
                row_colors.append(color_map.get(nucleotide, '#CCCCCC'))
                row_text.append(nucleotide)
            z_values.append(row_z)
            colors.append(row_colors)
            text_values.append(row_text)

        # Create the main heatmap with discrete colors
        fig_heatmap = go.Heatmap(
            z=z_values,
            x=positions,
            y=seq_names,
            text=text_values,
            texttemplate="%{text}",
            colorscale=self._create_discrete_colorscale(color_map),
            showscale=False,  # No legend
            hovertemplate="Sequence: %{y}<br>Position: %{x}<br>Nucleotide: %{text}<extra></extra>"
        )

        # Add consensus bar if requested
        if getattr(config, 'show_consensus', False):
            from .consensus import ConsensusCalculator

            # Use optimized sequences for consensus calculation to match the heatmap
            optimized_sequences = []
            for i, seq_name in enumerate(seq_names):
                # Reconstruct sequence from optimized dataframe
                seq_str = ''.join(df_optimized.iloc[i].tolist())
                # Create a temporary sequence object
                from .sequence import Sequence
                optimized_sequences.append(Sequence(seq_name, seq_str))

            consensus_calc = ConsensusCalculator(optimized_sequences)
            consensus_data = consensus_calc.calculate_alignment_consensus(
                ignore_gaps=getattr(config, 'consensus_ignore_gaps', True)
            )

            # Extract agreement percentages for bar chart
            agreement_values = [pos['agreement_percentage'] for pos in consensus_data]

            # Create consensus bar chart with no background
            fig_consensus = go.Bar(
                x=positions,
                y=agreement_values,
                marker_color='#666666',
                name='Consensus Agreement (%)',
                hovertemplate="Position: %{x}<br>Agreement: %{y:.1f}%<extra></extra>"
            )

            # Combine consensus bar and heatmap
            fig = make_subplots(
                rows=2, cols=1,
                shared_xaxes=True,
                row_heights=[0.25, 0.75],
                vertical_spacing=0.05,
            )

            fig.add_trace(fig_consensus, row=1, col=1)
            fig.add_trace(fig_heatmap, row=2, col=1)

            # Update layout for both subplots
            fig.update_layout(
                height=int(getattr(config, 'container_height', '400px').replace('px', '')),
                showlegend=False,
                plot_bgcolor='white',
                paper_bgcolor='white'
            )

            # Update axes with fixed range slider showing only 100 bases at a time
            # Calculate the initial visible range (first 100 bases or all if less than 100)
            visible_range = min(200, total_positions)
            initial_range = [0, visible_range - 1] if total_positions > 100 else [0, total_positions - 1]

            fig.update_yaxes(title_text="Agreement (%)", row=1, col=1)

            # Configure range slider for both consensus and heatmap
            fig.update_xaxes(
                rangeslider=dict(
                    visible=True,
                    thickness=0.1,
                    bgcolor='rgba(200, 200, 200, 0.3)',
                    bordercolor='rgba(100, 100, 100, 0.5)',
                    borderwidth=1
                ),
                range=initial_range,
                row=2, col=1
            )

            # Set the same range for the consensus subplot
            fig.update_xaxes(range=initial_range, row=1, col=1)

        else:
            # Single heatmap
            fig = go.Figure(data=[fig_heatmap])

            # Calculate the initial visible range (first 100 bases or all if less than 100)
            visible_range = min(200, total_positions)
            initial_range = [0, visible_range - 1] if total_positions > 100 else [0, total_positions - 1]

            fig.update_layout(
                xaxis_title='Position',
                yaxis_title='Sequence',
                xaxis=dict(
                    type='category',
                    rangeslider=dict(
                        visible=True,
                        thickness=0.1,
                        bgcolor='rgba(200, 200, 200, 0.3)',
                        bordercolor='rgba(100, 100, 100, 0.5)',
                        borderwidth=1
                    ),
                    range=initial_range
                ),
                yaxis=dict(type='category'),
                plot_bgcolor='white',
                paper_bgcolor='white',
                showlegend=False,
                margin=dict(l=40, r=40, t=40, b=40)
            )

            if hasattr(config, 'container_height'):
                fig.update_layout(height=int(config.container_height.replace('px', '')))

        return fig

    def _nucleotide_to_number(self, nucleotide):
        """Map nucleotides to numbers for heatmap z-values"""
        nuc_map = {'A': 0, 'G': 1, 'C': 2, 'T': 3, 'U': 3, '-': 4, 'N': 5}
        return nuc_map.get(nucleotide.upper(), 5)

    def _create_discrete_colorscale(self, color_map):
        """Create a discrete colorscale for nucleotides"""
        # Create colorscale based on nucleotide mapping
        colorscale = [
            [0.0, color_map.get('A', '#4c72a5')],   # A = 0
            [0.2, color_map.get('G', '#48a365')],   # G = 1
            [0.4, color_map.get('C', '#d0694a')],   # C = 2
            [0.6, color_map.get('T', '#e1c72f')],   # T = 3
            [0.8, color_map.get('-', '#f0f0f0')],   # Gap = 4
            [1.0, color_map.get('N', '#cccccc')]    # Unknown = 5
        ]

        return colorscale

    def _optimize_data_for_plotly(self, df, config: DisplayConfig):
        """Optimize data for Plotly rendering by applying decimation and chunking"""
        import warnings

        # Get original dimensions
        n_sequences, n_positions = df.shape

        # Check if optimization is needed and enabled
        needs_optimization = (
            n_positions > 5000  or
            n_sequences > 100
        )

        if not needs_optimization or not config.plotly_decimation:
            return df

        # Warn user about data decimation
        if needs_optimization:
            warnings.warn(
                f"Large alignment detected ({n_sequences} sequences × {n_positions} positions). "
                f"Applying data decimation for better performance. "
                f"Adjust ncols/nseqs to control this behavior.",
                UserWarning
            )

        # Apply sequence decimation if needed
        if n_sequences > config.plotly_max_sequences and config.plotly_max_sequences != 0:
            # Keep first and last sequences, then sample evenly from the middle
            keep_first = min(50, config.plotly_max_sequences // 4)
            keep_last = min(50, config.plotly_max_sequences // 4)
            keep_middle = config.plotly_max_sequences - keep_first - keep_last

            if keep_middle > 0:
                # Sample middle sequences
                middle_start = keep_first
                middle_end = n_sequences - keep_last
                middle_indices = list(range(middle_start, middle_end,
                                          max(1, (middle_end - middle_start) // keep_middle)))[:keep_middle]

                # Combine indices
                selected_indices = (list(range(keep_first)) +
                                  middle_indices +
                                  list(range(n_sequences - keep_last, n_sequences)))
            else:
                selected_indices = (list(range(keep_first)) +
                                  list(range(n_sequences - keep_last, n_sequences)))

            df = df.iloc[selected_indices]

        # Apply position decimation if needed
        if n_positions > config.plotly_max_positions and config.plotly_max_positions != 0:
            # Sample positions evenly across the alignment
            step_size = max(1, n_positions // config.plotly_max_positions)
            selected_positions = list(range(0, n_positions, step_size))[:config.plotly_max_positions]
            df = df.iloc[:, selected_positions]

        return df

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
                width: 1ch;
                background-color: grey;
                margin: 0;
                padding: 0;
                box-sizing: border-box;
                position: relative;
                cursor: pointer;
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
                consensus_html += '<div class="consensus-spacer"></div>'

            # Calculate bar height as percentage of container height
            height_percent = agreement  # agreement is already 0-100
            bar_html = (
                f'<div class="consensus-bar" '
                f'style="height: {height_percent}%;" '
                f'data-tooltip="Pos {pos}: {agreement:.1f}%">'
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
            ruler = self._create_ruler(padding, config.ncols, config.start_pos, config.block_size)
            print(ruler)

        for sequence in sequences:
            colored_seq = self.formatter.format_sequence(
                sequence.sequence,
                config.ncols,
                config.block_size,
                config.start_pos
            )
            print(f"{sequence.header:<{max_header_len}} {colored_seq}")

    def _create_ruler(self, padding: str, ncols: int, start_pos: int, block_size: int) -> str:
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

        padding_html = '&nbsp;' * len(padding)
        return padding_html + ''.join(numbers) + '\n' + padding_html + ticks
