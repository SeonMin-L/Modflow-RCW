import sys
rw_path = ""  # Input file paths for Radial_Well.py have been set
sys.path.append(rw_path)
import flopy
import Radial_Well as RW

sim_name = "" # Input modflow name
sim_path = "" # Input modflow folder path
exe_path = "" # Input mf6.exe path
pipe_excel = "" # Input pipe.excel path

sim = flopy.mf6.MFSimulation(sim_name=sim_name, version="mf6", exe_name="mf6.exe", sim_ws=sim_path)
iterative_data, Cassion_input_data, sim, RCW_out = RW.run(sim, pipe_excel, sim_name)


