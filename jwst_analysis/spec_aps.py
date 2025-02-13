import astropy.units as u
from photutils.aperture import EllipticalAperture, CircularAperture

"""Parameters for New and Combined MIRI/ NIRSPEC"""

# defined by ch1-long 
# parameters of left elliptical aperture
pos1 = (9.7, 16)
a1 = 4
b1 = 8.5
theta1 = (-5*u.degree).to(u.rad).value

# parameters of right elliptical aperture
pos2 = (18.5, 16.5)
a2 = 3.2
b2 = 8.3
theta2 = (-2*u.degree).to(u.rad).value

# parameters of circular aperture
pos = (14.1, 16.3)
r = 10.0

# create two elliptical apertures
e_ap1_new = EllipticalAperture(pos1, a1, b1, theta1)
e_ap2_new = EllipticalAperture(pos2, a2, b2, theta2)
# create circular aperutre
c_ap_new = CircularAperture(pos, r)

new_aps = [c_ap_new, e_ap1_new, e_ap2_new]



"""Parameters for old MIRI data"""

# defined by ch1-long
# parameters of left elliptical aperture
pos1 = (16.7, 20.2)
a1 = 4
b1 = 8.5
theta1 = (-5*u.degree).to(u.rad).value
# parameters of right elliptical aperture
pos2 = (24.8, 20.9)
a2 = 3.2
b2 = 8.3
theta2 = (-2*u.degree).to(u.rad).value

# parameters of circular aperture
pos = (20.5,20.5)
r = 10.0

#create elliptical apertures
e_ap1_old = EllipticalAperture(pos1, a1, b1, theta1)
e_ap2_old = EllipticalAperture(pos2, a2, b2, theta2)
# create circular aperutre
c_ap_old = CircularAperture(pos, r)

old_aps = [c_ap_old, e_ap1_old, e_ap2_old]




"""Parameters for JOYs MIRI Data"""
#defined in ch1-long
# parameters of left elliptical aperture
pos1 = (56.8, 20)
a1 = 4.2
b1 = 8.5
theta1 = (-5*u.degree).to(u.rad).value
# parameters of right elliptical aperture
pos2 = (65, 20)
a2 = 3.2
b2 = 8.3
theta2 = (-2*u.degree).to(u.rad).value
# parameters of circular aperture
pos=(60,20)
r =10.0

# create elliptical apertures
e_ap1_joys = EllipticalAperture(pos1, a1, b1, theta1)
e_ap2_joys = EllipticalAperture(pos2, a2, b2, theta2)
c_ap_joys = CircularAperture(pos, r)

joys_aps = [c_ap_joys, e_ap1_joys, e_ap2_joys]