import os

dir1 = '/Users/srikanthpatala/Desktop/cluster_structures/';
import pum_util_funcs as puf;

for file in os.listdir(dir1):
    if file.endswith(".out"):
        print(os.path.join(dir1, file))
        puf.analyze_gb_atoms(dir1, file)


