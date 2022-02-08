# Genetic markers that allow us to determine recombination break points
class GeneticMarker:
    def __init__(self, chrom, cm_map_pos, bp_map_pos)
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
    def __init__(self, pop_num, end_mark):
        self.pop_num = pop_num
        self.end_mark = end_mark
