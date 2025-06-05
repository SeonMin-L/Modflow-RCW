import sys
rw_path = ""  # Input file paths for Radial_Well.py and Gis_FloPy_nom.py have been set
sys.path.append(rw_path)
import Radial_Well as RW
import Gis_FloPy_nom as gf

model_setting_gis = "" # Input file path for the model_setting_gis
Modflow_work_file = "" # Input file path for modflow workstation
sim_name = "Test" # Input modflow name

sim = gf.modflow6_build(model_setting_gis, sim_name, Modflow_work_file)
pipe_excel= "" # Input pipe xlsx file path
iterative_data, Cassion_input_data, sim = RW.run(sim, pipe_excel, sim_name)

