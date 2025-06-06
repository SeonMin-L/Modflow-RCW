import flopy 
import os
import geopandas as gpd
import pickle
import numpy as np

def modflow6_build(Modflow_cell_path, sim_name, Modflow_work_file):
    
    with open(os.path.join(Modflow_cell_path, "MD_IP.txt"), "r") as file:
        lines = [line.rstrip("\n") for line in file.readlines()]
        nlay = int(lines[0].split(":")[1])
        nrow = int(lines[1].split(":")[1])
        ncol = int(lines[2].split(":")[1])
    
    chd_data=False
    riv_data=False
    drn_data=False
    
    kh = np.ones((nlay, nrow, ncol))
    kz = np.ones((nlay, nrow, ncol))
    Ss = np.ones((nlay, nrow, ncol))
    Sy = np.ones((nlay, nrow, ncol))
    InActive = np.ones((nlay, nrow, ncol))
    botm = np.ones((nlay, nrow, ncol))
    rch = np.zeros((nrow, ncol))
    riv = []
    drn = []
    chd = []
    for layer in range(nlay):
        Modflow_cell_load = gpd.read_file(os.path.join(Modflow_cell_path, f"MD_IP_Layer_{layer}.shp"))
        if layer == 0:
            top = np.array(Modflow_cell_load["up_el"]).reshape(nrow,ncol)        
            if "rch" in Modflow_cell_load.columns:
                rch = np.array(Modflow_cell_load["rch"]).reshape(nrow,ncol)
        kh[layer] = np.array(Modflow_cell_load["kh"]).reshape(nrow,ncol)
        kz[layer] = np.array(Modflow_cell_load["kv"]).reshape(nrow,ncol)
        Ss[layer] = np.array(Modflow_cell_load["ss"]).reshape(nrow,ncol)
        Sy[layer] = np.array(Modflow_cell_load["sy"]).reshape(nrow,ncol)
        InActive[layer] = np.array(Modflow_cell_load["inactive"]).reshape(nrow,ncol)
        botm[layer] = np.array(Modflow_cell_load["dn_el"]).reshape(nrow,ncol)
        
        if "riv" in Modflow_cell_load.columns:
            riv_data = True
            riv_in = Modflow_cell_load[Modflow_cell_load["riv"].notnull()]
            riv_in = [[int(data["layer"]), int(data["row"]), int(data["col"]), float(data["riv"].split(",")[0]), 
                    float(data["riv"].split(",")[1]), float(data["riv"].split(",")[2])] for idx, data in riv_in.iterrows()]
            riv = riv+riv_in
            
        if "drn" in Modflow_cell_load.columns:
            riv_data = True
            drn_in = Modflow_cell_load[Modflow_cell_load["drn"].notnull()]
            drn_in = [[int(data["layer"]), int(data["row"]), int(data["col"]), float(data["drn"].split(",")[0]), 
                    float(data["drn"].split(",")[1])] for idx, data in drn_in.iterrows()]
            drn = drn+drn_in
            
        if "chd" in Modflow_cell_load.columns:
            chd_data = True
            chd_in = Modflow_cell_load[Modflow_cell_load["chd"].notnull()]
            chd_in = [[int(data["layer"]), int(data["row"]), int(data["col"]), float(data["chd"].split(",")[0])] for idx, data in chd_in.iterrows()]
            chd = chd+chd_in
    
    row_array = Modflow_cell_load[(Modflow_cell_load["row"]==0)]["geometry"].bounds
    delc = np.array(row_array["maxx"]-row_array["minx"])
    col_array = Modflow_cell_load[(Modflow_cell_load["col"]==0)]["geometry"].bounds
    delr = np.array(col_array["maxy"]-col_array["miny"])

                    
    sim = flopy.mf6.MFSimulation(sim_name=sim_name, version="mf6", exe_name=os.path.join(Modflow_work_file, "mf6.exe"), sim_ws=Modflow_work_file)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=sim_name, save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc,top=top, botm=botm, idomain=InActive)
    
    strt = np.repeat(top, nlay, axis=0)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=kh, k33=kz)
    flopy.mf6.ModflowGwfsto(gwf, sy = Sy, ss=Ss, steady_state={0: True}) # transient={0: False}
    
    flopy.mf6.ModflowGwfrcha(gwf, recharge={0 :rch})
    
    if chd_data == True:
        chd_spd = {0:chd}
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    if riv_data == True:
        riv_spd = {0:riv}
        flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd)
    if drn_data == True:
        drn_spd = {0:drn}
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=drn_spd)
    
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",   
        outer_dvclose=1e-4,   
        outer_maximum=500,
        inner_maximum=100,
        inner_dvclose=1e-4,            
        rcloserecord=1e-5,             
        linear_acceleration="BICGSTAB",
        relaxation_factor=0.97         
    )
    
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f"{sim_name}.hds",budget_filerecord=f"{sim_name}.cbc",
                                saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])
    
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    
    return sim
