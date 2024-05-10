import numpy as np

root = "/home/alankar/Documents/cloud-crushing_PLUTO/output-wb-chi100.0eta100.0mach1.50tcoolmBtcc1.00e-02Tcl4.00e+04met1.00-w_cool-res8/clustering"
file_no = 20
data = np.load( f"{root}/cluster-labels.{file_no:04d}.npy" )

output.CellData.append(data[3], "label")
output.CellData.append(data[4], "probability")

