# Genetic markers that allow us to determine recombination break points
class GeneticMarker:
    def __init__(self, chrom, cm_map_pos, bp_map_pos):
        self.chrom = chrom
        self.cm_map_pos = cm_map_pos
        self.bp_map_pos = bp_map_pos

    def get_chrom(self):
        return self.chrom

    def get_map_pos(self):
        return self.cm_map_pos

    def get_bp_pos(self):
        return self.bp_map_pos


# stored segment lengths that allow us to recover population at marker locations
class HaplotypeSegment:
    def __init__(self, pop_num, chrom, end_coord):
        """
        Note the beginning of the segment is inferred based on previous
        segments stored throughout the simulation process. 
        Arguments
            pop_num - population label
            chrom - chromosome the haplotype segment lies on.
            end_coord - Ending coordinate of the haplotype segment
        """
        self.pop_num = pop_num
        self.chrom = chrom
        self.end_coord = end_coord

    def get_end_coord(self):
        return self.end_coord

    def get_chrom(self):
        return self.chrom

    def get_pop(self):
        return self.pop_num
