import flopy
import pandas as pd
import math
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

non_linear=True
convergence_loss=True
slot_loss=True
friction_loss= True
momemtom_loss= True
cassion_loss = True
Damping_start=True
caisson_inactive= True

u=1.0020 * 10**(-3)   #Ns/m2 = kg/m*s
p=1000   #kg/m
C=0.135
g=9.8*86400**2 #[m/d^2]
caisson_c = 1
convergence_end = 1
convergence_end_ratio = 0.01
cdh_end=0.05
iter_max = 100


# ✅ 필요 함수----------------------------------------------------------------------------------------------------
class Radial_well_input:
    def __init__(self, pipe_excel, gwf, sim):
        with pd.ExcelFile(pipe_excel) as excel_file:
            sheet_names = excel_file.sheet_names
            Cassion_input_data = {}
            for Cassion in sheet_names:
                input_CS={}
                Cassion_execl = pd.read_excel(excel_file, sheet_name=Cassion)
                input_CS["head"] = Cassion_execl.loc[0,"cassion"]
                input_CS["row"]= Cassion_execl.loc[1,"cassion"]
                input_CS["col"] = Cassion_execl.loc[2,"cassion"]
                input_CS["well_D"] = Cassion_execl.loc[3,"cassion"]
                input_CS["Skin_K"] = Cassion_execl.loc[4,"cassion"]
                input_CS["Skin_B"] = 0.5*(input_CS["Skin_K"]*100/86400)**(-3/4)*(100/86400)**2  #[day^2/m^2]
                input_CS["skin_Thick"] = Cassion_execl.loc[5,"cassion"]
                input_CS["pipe_open"] = Cassion_execl.loc[6,"cassion"]
                input_CS["ab_roughness"] = (Cassion_execl.loc[7,"cassion"]*10**(-3))/Cassion_execl.loc[3,"cassion"]
                input_CS["B_coefficient"] = Cassion_execl.loc[8,"cassion"]*0.001
                
                pipe_count = Cassion_execl.columns.str.contains('layer').sum()
                input_pipe = {}
                for pipe in range(pipe_count):
                    input_pipe[pipe] = Cassion_execl[[f"layer{pipe+1}", f"row{pipe+1}", f"col{pipe+1}", f"length{pipe+1}"]].dropna()
                    input_pipe[pipe].columns = ["layer", "row", "col", "length"]
                    input_pipe[pipe]["layer"] = input_pipe[pipe]["layer"].astype(int)
                    input_pipe[pipe]["row"] = input_pipe[pipe]["row"].astype(int)
                    input_pipe[pipe]["col"] = input_pipe[pipe]["col"].astype(int)
                    
                Cassion_input_data[Cassion] = [input_CS, input_pipe]
        
        self.gwf = gwf
        self.sim = sim
        self.Cassion_input_data = Cassion_input_data
        self.npf = gwf.get_package("npf")
        self.dis = gwf.get_package("dis")
        success, buff = sim.run_simulation(silent=True)
        headfile = flopy.utils.HeadFile(f"{gwf.model_ws}/{gwf.name}.hds")
        self.aqu_head = headfile.get_data(totim=1.0)
    
    def initial_drn_input(self, trn_hmd):
        idomain = self.dis.idomain.array
        initial_drn = []
        for Cassion in self.Cassion_input_data.keys():
            cassion_detail_data = self.Cassion_input_data[Cassion][0]
            cassion_head = cassion_detail_data["head"]
            rw = cassion_detail_data["well_D"]/2
            sk_thick = cassion_detail_data["skin_Thick"]
            B = cassion_detail_data["B_coefficient"]
            φ = cassion_detail_data["pipe_open"]
            rs = rw+sk_thick
            self.Cassion_input_data[Cassion][0]["rs"] = rs
            Ks = cassion_detail_data["Skin_K"]
            cassion_row = cassion_detail_data["row"]
            cassion_col = cassion_detail_data["col"]
            for key, Pipe in self.Cassion_input_data[Cassion][1].items():
                for idx, cell in Pipe.iterrows():
                    layer = int(cell["layer"] )
                    row = int(cell["row"])
                    col = int(cell["col"])
                    if trn_hmd:
                        head = (self.dis.top.array[row,col] + cassion_head)/2
                    else:
                        head = cassion_head
                
                    Kx = self.npf.k.array[layer, row, col]
                    Kz = self.npf.k33.array[layer, row, col]
                    dx = self.dis.delr.array[col]
                    botm = self.dis.botm.array[layer,row,col]
                    if layer == 0:
                        top = self.dis.top.array[row,col]
                    else:
                        top = self.dis.botm.array[layer-1,row,col]
                    dz = top-botm
                    re=0.28*math.sqrt((dx**2)*math.sqrt((Kz)/(Kx))+(dz**2)*math.sqrt((Kx)/(Kz)))/((Kz/Kx)**(1/4)+((Kx/Kz)**(1/4)))
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"re"] = re
                    Ka = self.Cassion_input_data[Cassion][1][key].loc[idx,"K_gm"] = math.sqrt(Kx*Kz)
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"B_gm"] = 0.5*(Ka*100/86400)**(-3/4)*(100/86400)**2  #[day^2/m^2]
                    Ks=Ka if sk_thick==0 else Ks
                    pipe_L = self.Cassion_input_data[Cassion][1][key].loc[idx,"length"]
                    CF = (1 / (2 * math.pi * pipe_L))
                    
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"S_aq"] = S_aq = CF * (1/(Ka)) * math.log(re/rs)
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"S_sk"] = S_sk = CF * (1/(Ks)) * math.log(rs/rw)
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"S_cv"] = S_cv = CF * (B/(rw*Ks*φ))
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"D_bt"] = (top+botm)/2
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"H"]    = head
                    
                    Fe = (S_aq+S_sk+S_cv) if convergence_loss == True else (S_aq+S_sk)
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"linear_CWC"] = linear_CWC = 1/Fe

                    drn_input = (layer, row, col, head, linear_CWC)
                    initial_drn.append(drn_input)
            
            if caisson_inactive == True:
                for layer_i in range(idomain.shape[0]):
                    idomain[int(layer_i), int(cassion_row), int(cassion_col)] = 0
        
        origin_drn=[]
        if "drn_0" in self.gwf.package_names:
            initial_drn_key = [(value[0], value[1], value[2]) for value in initial_drn]
            drn = self.gwf.get_package('drn_0')
            drn_data = drn.stress_period_data.data.get_data(0)
            for i in drn_data[0]:
                key = (i["cellid"][0],i["cellid"][1],i["cellid"][2])
                if key not in initial_drn_key:
                    origin_drn.append([i["cellid"][0],i["cellid"][1],i["cellid"][2],i["elev"],i["cond"]])
        
        initial_drn = origin_drn+initial_drn
        self.gwf.dis.idomain.set_data(idomain)
        if "drn_0" in self.gwf.package_names:
            self.gwf.remove_package("drn_0")
        flopy.mf6.ModflowGwfdrn(self.gwf, stress_period_data={0: initial_drn}, save_flows=True)
        self.sim.write_simulation(silent=True)
        
        return self.gwf, origin_drn

def output_read(Cassion_input_data, repeat, gwf):
    Q_in_max =0
    nlay = gwf.dis.nlay.get_data()
    nrow = gwf.dis.nrow.get_data()
    ncol = gwf.dis.ncol.get_data()
    def node_equation(node, nlay, nrow, ncol):
        layer = node//(nrow*ncol)
        remain_cell = node % (nrow*ncol)
        row = remain_cell//ncol
        col = remain_cell % ncol
        return layer, row, col
    
    cbc = flopy.utils.CellBudgetFile(f"{gwf.model_ws}/{gwf.name}.cbc", precision="double")
    drain_leakage = cbc.get_data(text="DRN")
    Drn_out_df = pd.DataFrame({
        'layer': [node_equation(node-1, nlay, nrow, ncol)[0] for node in drain_leakage[0]["node"]],
        'row': [node_equation(node-1, nlay, nrow, ncol)[1] for node in drain_leakage[0]["node"]],
        'col': [node_equation(node-1, nlay, nrow, ncol)[2] for node in drain_leakage[0]["node"]],
        'leakage': drain_leakage[0]["q"]})
    
    for cassion_key, cassion_data in Cassion_input_data.items():
        for pipe_key, pipe_data in cassion_data[1].items():
            Q_sum = 0
            for idx, cell in pipe_data.iloc[::-1].iterrows():
                row = int(cell["row"])
                col = int(cell["col"])
                layer = int(cell["layer"])
                target_row = Drn_out_df[(Drn_out_df["layer"]==layer) & (Drn_out_df["row"]==row) & (Drn_out_df["col"]==col)].index
                leakage_value = Drn_out_df.loc[target_row, "leakage"].values[0]
                Q_sum += leakage_value*-1
                
                if repeat >= 1:
                    old_Q =  Cassion_input_data[cassion_key][1][pipe_key].loc[idx, "Q_in"]
                    new_Q = leakage_value*-1
                    Q_change = abs(old_Q - new_Q)
                    Q_in_max = max(Q_in_max, Q_change)
                Cassion_input_data[cassion_key][1][pipe_key].loc[idx, "Q_in"] = leakage_value*-1
                Cassion_input_data[cassion_key][1][pipe_key].loc[idx, "Q_pipe"] = Q_sum
    
    headfile = flopy.utils.HeadFile(f"{gwf.model_ws}/{gwf.name}.hds")
    aqu_head = headfile.get_data(totim=1.0)
    
    return Cassion_input_data, aqu_head, Q_in_max

class Horizontal_calculate:
    def __init__(self, Cassion_data, aqu_head):
        self.head = Cassion_data[0]["head"]
        self.Skin_K = Cassion_data[0]["Skin_K"]
        self.Skin_B = Cassion_data[0]["Skin_B"]
        self.rw = Cassion_data[0]["well_D"]/2
        self.rs = Cassion_data[0]["rs"]
        self.ab_roughness = Cassion_data[0]["ab_roughness"]
        self.pipe_data_all = Cassion_data[1]
        self.skin_Thick = Cassion_data[0]["skin_Thick"]
        self.pipe_open = Cassion_data[0]["pipe_open"]
        self.aqu_head = aqu_head
        
    def pipe(self, top, alpha):
        def f_calculation(e, pV, D, iq, p, ReT=3000, ReL=2300):     # Colebrook-White equation
            pV = np.asarray(pV, dtype=float)
            iq = np.asarray(iq, dtype=float)
            Re=p*pV*D/u                          # [kg/m^3]*[m/s]*[m]/[kg/m*s] = []
            WRe = p*iq/(math.pi*u)               # [kg/m^3]*[m*2/s]/[kg/m*s] = []
            f = np.zeros_like(pV, dtype=float)
            nz = (pV != 0) & np.isfinite(Re) & np.isfinite(WRe)
            
            m_turb = nz & (Re >= ReT)
            m_lam  = nz & (Re <= ReL)
            m_mid  = nz & (~m_turb) & (~m_lam)  # ReL < Re < ReT
            
            def f_lam(Re, WRe):
                return 64/Re*(1+0.04304*WRe**0.6142)
                
            def f_turb(Re, WRe):
                f1=(1/(-2 * np.log10(e/3.7065 - (5.0272/Re) * np.log10(e/3.827 - 
                   (4.567/Re)*np.log10((e/7.7918)**0.9924 + (5.3326/(208.5 + Re))**0.9345)))))**2
                return f1*(1-0.0153*WRe**0.3978)
            
            if np.any(m_turb):
                f[m_turb] = f_turb(Re[m_turb], WRe[m_turb])
            if np.any(m_lam):
                f[m_lam] = f_lam(Re[m_lam], WRe[m_lam])    
            if np.any(m_mid):
                x2 = f_turb(np.full(np.count_nonzero(m_mid), ReT, dtype=float), WRe[m_mid])
                x1 = f_lam (np.full(np.count_nonzero(m_mid), ReL, dtype=float), WRe[m_mid])
                f[m_mid] = (x2 - x1) / (ReT - ReL) * (Re[m_mid] - ReL) + x1
            
            return f
        
        def f_head_L(IQ, SQ, f2, f1, L, D):
            f = 0.5 * (f1 + f2)
            Q2_int = L * ( (IQ**2)/3.0 + IQ*SQ + SQ**2 )
            Cf = 8.0 / (np.pi**2 * g * D**5)
            dhf = Cf * f * Q2_int
            return dhf
        
        def m_head_L(IQ, SQ, L, D):
            Cm = 8.0 / (np.pi**2 * g * D**4)
            dhm = Cm * (IQ**2 + 2*IQ*SQ)
            return dhm
    
        def c_head_L(V, f1):
            dhc = np.zeros_like(f1)
            dhc[0] = V**2 / (2*9.8)
            return dhc
            
        mod_input  = []
        rel_input = []
        for pipe_key, pipe_data in self.pipe_data_all.items():
            L  = self.pipe_data_all[pipe_key]["length"].to_numpy()
            pipe_data["Q_inp"]  = (pipe_data["Q_in"] / L)/86400  #m2/sec
            pipe_data["V_pipe"] = (pipe_data["Q_pipe"] / (math.pi * self.rw ** 2))/86400              #([m^3/d]/[m^2])*(d/s)=m/s
            pipe_data["f"]      = f_calculation(self.ab_roughness, pipe_data["V_pipe"], self.rw*2, pipe_data["Q_inp"], p)
            
            SQ = np.concatenate([pipe_data["Q_pipe"][1:].to_numpy(), [0.0]])
            IQ = pipe_data["Q_in"].to_numpy()
            f2 = pipe_data["f"].to_numpy()
            f1 = np.concatenate([pipe_data["f"][1:].to_numpy(), [0.0]])
            CF = (1 / (2 * math.pi * L)) ** 2
            
            pipe_data["F_Head_L"] = Fa    = 0 if not friction_loss else f_head_L(IQ, SQ, f2, f1, L, self.rw*2)
            pipe_data["M_Head_L"] = Ma    = 0 if not momemtom_loss else m_head_L(IQ, SQ, L, self.rw*2)
            pipe_data["C_Head_L"] = Ca    = 0 if not cassion_loss  else c_head_L(pipe_data.loc[0,"V_pipe"], f1)
            pipe_data["S_aq_N"]   = NS_aq = 0 if not non_linear    else CF * pipe_data["B_gm"] * (1 / self.rs - 1 / pipe_data["re"]) * pipe_data["Q_in"]
            pipe_data["S_sk_N"]   = NS_sk = 0 if not non_linear    else CF * self.Skin_B * (1 / self.rw - 1 / self.rs) * pipe_data["Q_in"]
            pipe_data["S_sl_N"]   = NS_sl = 0 if not slot_loss     else CF * (1 / (2 * g)) * (1 / (self.rw * C * self.pipe_open)) ** 2 * pipe_data["Q_in"]

            for i, __ in enumerate(Fa):
                l, r, c = pipe_data.loc[i,"layer"], pipe_data.loc[i,"row"], pipe_data.loc[i,"col"]
                ah, th  = self.aqu_head[l,r,c], top[r,c]
                NH   = self.head + Ca[0] + (sum(Fa[0:i]) + sum(Ma[0:i])) + (Fa[i] + Ma[i])/2
                OH   = pipe_data.loc[i,"H"]
                CWC  = 1 / (1/pipe_data.loc[i,"linear_CWC"] + NS_aq[i] + NS_sk[i] + NS_sl[i])
                
                if np.isnan(NH) or NH > th:
                    NH = th
                aeff = alpha
                UH = (1-aeff) * OH + aeff * NH
                
                pipe_data.loc[i, "H"]    = UH
                pipe_data.loc[i, "CWC"]  = CWC
                pipe_data.loc[i, "RH"]   = NH
            
            mod_input.extend([[int(i["layer"]), int(i["row"]), int(i["col"]), i["H"],  i["CWC"]] for __, i in pipe_data.iterrows()])
            rel_input.extend([[int(i["layer"]), int(i["row"]), int(i["col"]), i["RH"], i["CWC"]] for __, i in pipe_data.iterrows()])
            self.pipe_data_all[pipe_key] = pipe_data
            
        return self.pipe_data_all, mod_input, rel_input

def iterative_data_save(Cassion_input_data, iterative_data, repeat):
    iterative_data[repeat] = {}
    for Cassion_key, Cassion_data in Cassion_input_data.items():
        iterative_data[repeat][Cassion_key] = {}
        for pipe_key, pipe_df in Cassion_data[1].items():
            iterative_data[repeat][Cassion_key][pipe_key] = pipe_df.copy()
    return iterative_data
        
def model_run(sim, gwf, drn_input):
    gwf.remove_package("drn_0")
    flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: drn_input}, save_flows=True)
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    if not success:
        print("MODFLOW did not terminate normally.")
        
# ✅ 실행---------------------------------------------------------------------------------------------------------
def run(sim, pipe_excel, sim_name, alpha=0.3, cdh_end=0.5, cdq_end1 = 0.1, cdq_end2=5, iter_max = 100, trn_hmd = True):

    gwf = sim.get_model(sim_name)
    iterative_data = {}
    repeat = 0

    root = Radial_well_input(pipe_excel, gwf, sim)
    gwf, origin_drn = root.initial_drn_input(trn_hmd)
    success, buff = sim.run_simulation(silent=True)
    if not success:
        print("MODFLOW did not terminate normally.")
    Cassion_input_data, aqu_head, Q_in_max = output_read(root.Cassion_input_data, repeat, gwf)
    iterative_data = iterative_data_save(Cassion_input_data, iterative_data, repeat)
    top = gwf.dis.top.array
    
    while trn_hmd:
        repeat +=1
        total_input = []
        final_input = []
        for cassion_name, Cassion_data in Cassion_input_data.items():
            root = Horizontal_calculate(Cassion_data, aqu_head)
            Cassion_input_data[cassion_name][1], mod_input, rel_input = root.pipe(top, alpha)
            total_input.extend(mod_input)
            final_input.extend(rel_input)
            
        iterative_data = iterative_data_save(Cassion_input_data, iterative_data, repeat)
        drn_input = [[int(r[0]), int(r[1]), int(r[2]), *r[3:]] for r in (origin_drn + total_input)]
        model_run(sim, gwf, drn_input)
        Cassion_input_data, aqu_head, Q_in_max = output_read(Cassion_input_data, repeat, gwf)
        
        if repeat == 1:
            f_Q_in_max = Q_in_max
            print(f"{repeat} {Q_in_max: .3f}")
        elif repeat >= 2:
            now_input = np.array(total_input)
            rel_input = np.array(final_input)
            dh = max(abs(now_input[:,3] - rel_input[:,3]))
            dq = Q_in_max / f_Q_in_max
            print(f"{repeat} {dq: .3f}, {dh : .3f}, {Q_in_max : .3f}")
            if (cdh_end > dh and cdq_end1 > dq)  or (Q_in_max < cdq_end2):
                drn_input = [[int(r[0]), int(r[1]), int(r[2]), *r[3:]] for r in (origin_drn + final_input)]
                model_run(sim, gwf, drn_input)
                Cassion_input_data, aqu_head, Q_in_max = output_read(Cassion_input_data, repeat+1, gwf)
                break
            if repeat == iter_max:
                drn_input = []
                for cn, cdata in Cassion_input_data.items():
                    for pn, pdata in cdata[1].items():
                        H = np.zeros([pdata.shape[0],5])
                        H[:,0] = pdata["layer"].to_numpy()
                        H[:,1] = pdata["row"].to_numpy()
                        H[:,2] = pdata["col"].to_numpy()
                        for i in range(int(iter_max*0.9), iter_max):
                            H[:,3] += iterative_data[i][cn][pn].loc[:,"H"].to_numpy()
                            H[:,4] += iterative_data[i][cn][pn].loc[:,"CWC"].to_numpy()
                        Cassion_input_data[cn][1][pn]["H"]   = H[:,3] = H[:,3] / (iter_max - int(iter_max*0.9))
                        Cassion_input_data[cn][1][pn]["CWC"] = H[:,4] = H[:,4] / (iter_max - int(iter_max*0.9))
                        pipe_input = [[int(r[0]), int(r[1]), int(r[2]), *r[3:]] for r in H]
                        drn_input.extend(pipe_input)
                drn_input.extend(origin_drn)
                model_run(sim, gwf, drn_input)
                Cassion_input_data, aqu_head, Q_in_max = output_read(Cassion_input_data, repeat+1, gwf)
                break
    
    gwf.remove_package("drn_0")
    if len(origin_drn) > 0:
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: origin_drn}, save_flows=True)
    sim.write_simulation(silent=True)
    
    RCW_out={}
    for caisson, data in Cassion_input_data.items():
        Q_out = 0
        for idx, lateral in data[1].items():
            Q = lateral.loc[0,"Q_pipe"]
            Q_out += Q
        RCW_out[caisson] = Q_out
            
    return iterative_data, Cassion_input_data, sim, RCW_out, aqu_head
        
        
