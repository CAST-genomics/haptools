from karyogram import plot_karyogram

def main():
    sample_file = "/storage/mlamkin/data/simwas/ten_gen.bp" 
    title = "Ten Generations"
    centromeres = []
    out = "/storage/mlamkin/data/simwas/ten_gen.png"
    plot_karyogram(sample_file, title, centromeres, out)
    return

if __name__ == '__main__':
    main()
