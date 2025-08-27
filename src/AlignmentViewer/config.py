from typing import Optional

class DisplayConfig:
    """Configuration for alignment display"""
    def __init__(
        self,
        nseqs: int = 0,
        ncols: int = 0,
        show_ruler: bool = True,
        show_consensus: bool = False,
        consensus_height: str = "50px",
        consensus_ignore_gaps: bool = True,
        block_size: int = 10,
        start_pos: int = 0,
        container_height: str = "300px",
        as_html: Optional[bool] = None,
        as_plotly: Optional[bool] = None,
        plotly_decimation: bool = True,    # Enable automatic data decimation for large datasets
        plotly_chunk_size: int = 1000,     # Size of chunks for data processing,
        color_snps_only: bool = False
    ):
        self.nseqs = int(nseqs)
        self.ncols = int(ncols)
        self.show_ruler = show_ruler
        self.show_consensus = show_consensus
        self.consensus_height = consensus_height
        self.consensus_ignore_gaps = consensus_ignore_gaps
        self.block_size = int(block_size)
        self.start_pos = int(start_pos)
        self.container_height = container_height
        self.as_html = as_html
        self.as_plotly = as_plotly
        # Performance optimization parameters
        self.plotly_max_positions = int(ncols)
        self.plotly_max_sequences = int(nseqs)
        self.plotly_decimation = plotly_decimation
        self.plotly_chunk_size = int(plotly_chunk_size)
        self.color_snps_only = color_snps_only

    def validate(self):
        """Validate the configuration values."""
        if self.nseqs < 0:
            raise ValueError("Number of sequences must be greater than or equal to 0.")
        if self.ncols < 0:
            raise ValueError("Number of columns must be greater than or equal to 0.")
        if self.block_size < 3:
            raise ValueError("Block size must be greater than or equal to 3.")
        if self.start_pos < 0:
            raise ValueError("Start position must be greater than or equal to 0.")
        if self.plotly_chunk_size < 100:
            raise ValueError("plotly_chunk_size must be at least 100.")

        # Warn if both as_html and as_plotly are set to True
        if self.as_html is True and self.as_plotly is True:
            import warnings
            warnings.warn("Both as_html and as_plotly are set to True. as_plotly will take precedence.", UserWarning)
