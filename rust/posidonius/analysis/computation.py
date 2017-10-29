import numpy as np

#-------------------------------------------------------------------------------
# Main process
#-------------------------------------------------------------------------------
# The 2 eccentricity dependant factors in the equation in a
def _Na1(e):
    return (1.0+31.0/2.0*(e*e)+255.0/8.0*np.power(e, 4) \
            + 185.0/16.0* np.power(e, 6) + 25.0/64.0* np.power(e,8)) / np.power((1.0- (e*e)), (15.0/2.0))

def _Na2(e):
    return (1.+15./2.0*(e*e) + 45./8.0 * np.power(e, 4) + 5./16.0* np.power(e, 6)) / np.power((1.0- (e*e)), 6)

def _No2(e):
    return (1.0 + 3.0 * (e*e) + 3.0/8.0 * np.power(e, 4)) / np.power((1.0- (e*e)), 5)

# Mean orbital angular velocity without the a dependance         (m^3/2.s-1)
def _norb(G, Mp, Ms):
    return np.sqrt(G) * np.sqrt(Mp+Ms)

# Jeremys Ki factor :
def _Kplan(k2deltat_plan, G, Mp, Ms, Rp, a):
    return 3./2. * k2deltat_plan * (G*(Mp*Mp)/Rp) * np.power(Ms/Mp,2) \
            * np.power(Rp/a, 6) * np.power(_norb(G,Mp,Ms) * np.power(a, -1.5), 2)

def energydot(a, e, rotp, oblp, G, Mp, Ms, Rp, k2deltat_plan):
    return 2. * _Kplan(k2deltat_plan, G, Mp, Ms, Rp, a) \
            * (_Na1(e) - 2.0*_Na2(e) * np.cos(oblp) * (rotp/(_norb(G, Mp, Ms)*np.power(a, -1.5))) \
            + (1.0 + np.power(np.cos(oblp),2))/2.0 * _No2(e) * np.sqrt(1.0-(e*e)) * np.power(rotp/(_norb(G, Mp, Ms)* np.power(a, -1.5)), 2))

