# openCFS-2D-geometry

Optimizer used: ipopt,
method: SIMP,
Boundary Condition at Nozzle_curve: Pressure equation,

running Terminal Commands:
1) cfs -m box2d-t_0.3-nx_200-ny_20-nz_20.mesh 2D_distributed_load
2) show_density 2D_distributed_load.density.xml
3) plotviz.py 2D_distributed_load.plot.dat -x 1 -y 2 4
4) For comparison between multiple .dat files: plotviz.py *.plot.dat -x 1 -y 2 4
