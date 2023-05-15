import swiftest
import numpy as np

swiftersim = swiftest.Simulation(simdir="swifter_sim", read_data=True, codename="Swifter")
swiftestsim = swiftest.Simulation(simdir="swiftest_sim",read_data=True)
swiftestsim.data = swiftestsim.data.swap_dims({"name" : "id"})
swiftdiff = swiftestsim.data - swiftersim.data
swiftdiff['rh'].plot(x="time",hue="id",col="space")
swiftdiff['vh'].plot(x="time",hue="id",col="space")