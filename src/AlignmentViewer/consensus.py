from typing import List, Dict, Tuple
from collections import Counter
from .sequence import Sequence


class ConsensusCalculator:
    """Calculate consensus and agreement statistics for sequence alignments"""
    
    def __init__(self, sequences: List[Sequence]):
        self.sequences = sequences
        self.alignment_length = len(sequences[0].sequence) if sequences else 0
        self._validate_sequences()
    
    def _validate_sequences(self) -> None:
        """Validate that all sequences have the same length"""
        if not self.sequences:
            raise ValueError("No sequences provided")
        
        lengths = [len(seq.sequence) for seq in self.sequences]
        if not all(length == lengths[0] for length in lengths):
            raise ValueError("All sequences must have the same length for alignment")
    
    def calculate_position_consensus(self, position: int, ignore_gaps: bool = True) -> Dict[str, float]:
        """
        Calculate consensus statistics for a specific position
        
        Args:
            position: The position in the alignment (0-based)
            ignore_gaps: Whether to ignore gap characters in consensus calculation
            
        Returns:
            Dictionary with consensus statistics including agreement percentage
        """
        if position < 0 or position >= self.alignment_length:
            raise ValueError(f"Position {position} out of range (0-{self.alignment_length-1})")
        
        # Get all characters at this position
        chars = [seq.sequence[position].upper() for seq in self.sequences]
        
        # Filter out gaps if requested
        if ignore_gaps:
            chars = [char for char in chars if char not in ['-', '.', 'X']]
        
        if not chars:
            return {
                'most_common': '-',
                'agreement_percentage': 0.0,
                'total_sequences': len(self.sequences),
                'non_gap_sequences': 0,
                'character_counts': {}
            }
        
        # Count character frequencies
        counter = Counter(chars)
        most_common_char, most_common_count = counter.most_common(1)[0]
        
        # Calculate agreement percentage
        agreement_percentage = (most_common_count / len(chars)) * 100
        
        return {
            'most_common': most_common_char,
            'agreement_percentage': agreement_percentage,
            'total_sequences': len(self.sequences),
            'non_gap_sequences': len(chars),
            'character_counts': dict(counter)
        }
    
    def calculate_alignment_consensus(self, ignore_gaps: bool = True) -> List[Dict[str, float]]:
        """
        Calculate consensus for all positions in the alignment
        
        Args:
            ignore_gaps: Whether to ignore gap characters in consensus calculation
            
        Returns:
            List of consensus dictionaries, one per position
        """
        return [
            self.calculate_position_consensus(pos, ignore_gaps)
            for pos in range(self.alignment_length)
        ]
    
    def get_consensus_agreement_percentages(self, ignore_gaps: bool = True) -> List[float]:
        """
        Get agreement percentages for all positions (useful for bar chart visualization)
        
        Args:
            ignore_gaps: Whether to ignore gap characters in consensus calculation
            
        Returns:
            List of agreement percentages (0-100) for each position
        """
        return [
            self.calculate_position_consensus(pos, ignore_gaps)['agreement_percentage']
            for pos in range(self.alignment_length)
        ]
    
    def get_consensus_sequence(self, ignore_gaps: bool = True) -> str:
        """
        Generate a consensus sequence using the most common character at each position
        
        Args:
            ignore_gaps: Whether to ignore gap characters in consensus calculation
            
        Returns:
            Consensus sequence string
        """
        consensus_chars = []
        for pos in range(self.alignment_length):
            consensus_data = self.calculate_position_consensus(pos, ignore_gaps)
            consensus_chars.append(consensus_data['most_common'])
        
        return ''.join(consensus_chars)
    
    def generate_consensus_bar_data(self, start_pos: int = 0, ncols: int = 0, 
                                   block_size: int = 10, ignore_gaps: bool = True) -> List[Tuple[int, float]]:
        """
        Generate data for consensus bar chart visualization
        
        Args:
            start_pos: Starting position for the view
            ncols: Number of columns to display (0 for all)
            block_size: Size of sequence blocks for alignment with visualization
            ignore_gaps: Whether to ignore gap characters
            
        Returns:
            List of (position, agreement_percentage) tuples
        """
        # Get agreement percentages for the specified range
        end_pos = start_pos + ncols if ncols > 0 else self.alignment_length
        end_pos = min(end_pos, self.alignment_length)
        
        bar_data = []
        for pos in range(start_pos, end_pos):
            agreement = self.calculate_position_consensus(pos, ignore_gaps)['agreement_percentage']
            bar_data.append((pos, agreement))
        
        return bar_data
