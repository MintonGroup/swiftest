"""
Calculate the approx. mass of a ring given the central body observed parameters

"""
import numpy as np

simdir = '/scratch/bell/anand43/swiftest_runs/'

# ring parameters (in SI units)
density = 2000 # kg / m^3
width = 3.5e3 # m
radius = 404.8e3 # m # assume inner radius
outer_rad = radius + width
height = 0.1e3 # assume height of 100 m thick

# central body parameters (in SI units)
cb_mass = 3.874e20 # kg
cb_name = 'Chariklo'

vol = np.pi * (outer_rad**2 - radius**2) * height # ring volume
ring_mass = vol * density # ring mass

print_msg = f'\nMass of Ring around {cb_name} = {ring_mass} kg = {ring_mass / cb_mass} M_{cb_name}'
print(print_msg)
print(f'Assuming a ring thickness of {height} m')

file = open(simdir + cb_name + f'/{cb_name}_ring_mass.txt', 'w')
file.write(print_msg)