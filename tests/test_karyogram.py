import pytest
from pathlib import Path

from haptools.karyogram import GetHaplotypeBlocks, PlotKaryogram

DATADIR = Path(__file__).parent.joinpath("data")

def test_GetHaplotypeBlocks():
	test_file = DATADIR.joinpath("test.bp")
	sample_blocks = GetHaplotypeBlocks(test_file, "Sample_1")
	assert(sample_blocks[0][0]["pop"]=="YRI")
	assert(sample_blocks[0][0]["chrom"]==1)
	assert(sample_blocks[0][0]["start"]==0.0001)
	assert(sample_blocks[0][0]["end"]==168.003442)

	assert(sample_blocks[1][0]["pop"]=="YRI")
	assert(sample_blocks[1][0]["chrom"]==1)
	assert(sample_blocks[1][0]["start"]==0.0001)
	assert(sample_blocks[1][0]["end"]==87.107755)
	
	sample_blocks = GetHaplotypeBlocks(test_file, "Sample_2")
	assert(sample_blocks[0][-1]["pop"]=="YRI")
	assert(sample_blocks[0][-1]["chrom"]==2)
	assert(sample_blocks[0][-1]["start"]==180.837755+0.0001)
	assert(sample_blocks[0][-1]["end"]==244.341689)

	assert(sample_blocks[1][0]["pop"]=="YRI")
	assert(sample_blocks[1][0]["chrom"]==1)
	assert(sample_blocks[1][0]["start"]==0.0001)
	assert(sample_blocks[1][0]["end"]==85.107755)

test_file = DATADIR.joinpath("5gen.bp")
centromeres = DATADIR.joinpath("centromeres_hg19.txt")
PlotKaryogram(test_file, "Sample_1", "test.png",
       centromeres_file=centromeres, title=None)