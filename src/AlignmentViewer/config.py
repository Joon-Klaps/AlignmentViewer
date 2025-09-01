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
