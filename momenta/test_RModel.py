

import numpy as num
from radial import RadialModel
from pylab import figure,plot,show,imshow,close
from pdb import set_trace as stop

radii = num.arange(100.)
intens = 100. * num.exp(-radii/20.)

center = (100.,100.)
dims = (200,300)
q = 0.5
pa = 15.

Model = RadialModel(radii,intens,dims,center,q,pa,c=0)

figure(2)
imshow(Model,origin='lower')
show()
close(2)


stop()

