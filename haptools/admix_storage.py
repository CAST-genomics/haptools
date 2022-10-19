# Genetic markers that allow us to determine recombination break points
class GeneticMarker:
    def __init__(self, chrom, cm_map_pos, bp_map_pos, prev_coord):
        self.chrom = chrom
        self.cm_map_pos = cm_map_pos
        self.bp_map_pos = bp_map_pos
        self.prev_coord = prev_coord

    def __repr__(self):
        return f"(Chrom {self.chrom}, {self.cm_map_pos} cM, {self.bp_map_pos} bp)"

    def __str__(self):
        return f"(Chrom {self.chrom}, {self.cm_map_pos} cM, {self.bp_map_pos} bp)"

    def get_chrom(self):
        return self.chrom

    def get_map_pos(self):
        return self.cm_map_pos

    def get_bp_pos(self):
        return self.bp_map_pos

    def get_prev_coord(self):
        return self.prev_coord


# stored segment lengths that allow us to recover population at marker locations
class HaplotypeSegment:
    def __init__(self, pop_num, chrom, end_coord, end_pos):
        """
        Note the beginning of the segment is inferred based on previous
        segments stored throughout the simulation process.

        Parameters
        ----------
            pop_num
                population label
            chrom
                chromosome the haplotype segment lies on.
            end_coord
                Ending coordinate in bp of the haplotype segment
            end_pos
                Ending coordinate in centimorgans of the hap segment
        """
        self.pop_num = pop_num
        self.chrom = chrom
        self.end_coord = end_coord
        self.end_pos = end_pos

    def __repr__(self):
        return (
            f"(Population {self.pop_num}, chrom {self.chrom}, "
            + f"end_coord {self.end_coord}, end_pos {self.end_pos})"
        )

    def __str__(self):
        return (
            f"Population {self.pop_num}, chrom {self.chrom}, "
            + f"end_coord {self.end_coord}, end_pos {self.end_pos}"
        )

    def get_end_pos(self):
        return self.end_pos

    def get_end_coord(self):
        return self.end_coord

    def get_chrom(self):
        return self.chrom

    def get_pop(self):
        return self.pop_num
