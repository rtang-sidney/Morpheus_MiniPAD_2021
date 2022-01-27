import numpy as np

COIL_O2 = 'tti1_out1'
COIL_I1 = 'tti1_out2'
COIL_I2 = 'tti2_out1'
COIL_O1 = 'tti2_out2'
CHANNELS = [COIL_I1, COIL_O1, COIL_I2, COIL_O2]

I1_POSITION = 'upstream, inner'
O1_POSITION = 'upstream, outer'
I2_POSITION = 'downstream, inner'
O2_POSITION = 'downstream, outer'
POSITIONS = [I1_POSITION, O1_POSITION, I2_POSITION, O2_POSITION]

COILS_POSITIONS = dict(zip(CHANNELS, POSITIONS))
