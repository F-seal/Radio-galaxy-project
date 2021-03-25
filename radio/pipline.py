from RadioData import *

for i in (5,10,15,20,30,40,50,75):
    rd = RadioData(flux_cut=i,catalog_name='nvss')
    rd.get_cl()
    rd.print_info()
    rd.get_distribution()
    rd.save_map(mask=True)
    rd.save_map(mask=True)
    rd.plot_cl()