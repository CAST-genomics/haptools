from karyogram import plot_karyogram

def main():
    sample_name = "Sample_1"
    sample_file = "/storage/mlamkin/data/simwas/ten_gen.bp" 
    title = "Ten Generations"
    chrX = False
    centromeres = []
    pop_order = ['1','2']
    colors = []
    out = "/storage/mlamkin/data/simwas/ten_gen.png"
    plot_karyogram(sample_name, sample_file, title, chrX, 
                   centromeres, pop_order, colors, out)
    return

if __name__ == '__main__':
    main()
