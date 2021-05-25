import math
import numpy as np

def rotatePC(ins,outs):

    x = ins['X']
    y = ins['Y']

    # Transform to FRF coords (state plan to frf) using the technique in frfCoord.m script #
    spAngle = (90-69.974707831)/(180/math.pi) # Angle of FRF system relative to state plane #
    Eom = 901951.6805 # FRF origin state plane Easting #
    Nom = 274093.1562 # FRF origin state plane Northing #

    SpLengE = x-Eom
    SpLengN = y-Nom
    R = np.sqrt(SpLengE**2+SpLengN**2)
    Ang1 = np.arctan2(SpLengE,SpLengN)
    Ang2 = Ang1+spAngle
    X_rot = np.multiply(R,np.sin(Ang2))
    Y_rot = np.multiply(R,np.cos(Ang2))

    outs['X'] = X_rot
    outs['Y'] = Y_rot

    return True
