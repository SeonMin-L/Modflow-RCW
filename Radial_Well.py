import flopy
import pandas as pd
import math
import sympy as sp

non_linear=True
convergence_loss = True
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
alpha=0.25
convergence_end = 1
convergence_end_ratio = 0.03

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
    
    def initial_drn_input(self):
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
                    
                    if sk_thick==0:
                        rs=rw
                        Ks=Ka
                    
                    pipe_L = self.Cassion_input_data[Cassion][1][key].loc[idx,"length"]
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"S_aq"] = S_aq =(1/(Ka)) * math.log(re/rs)
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"S_sk"] = S_sk=(1/(Ks)) * math.log(rs/rw)
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"S_cv"] = S_cv=(B/(rw*Ks*φ))
                    
                    if convergence_loss == True:               
                        Fe=(1/(2*math.pi*pipe_L))*(S_aq+S_sk+S_cv)
                    else:
                        Fe=(1/(2*math.pi*pipe_L))*(S_aq+S_sk)
                        
                    self.Cassion_input_data[Cassion][1][key].loc[idx,"linear_CWC"] = linear_CWC = 1/Fe
                    
                    drn_input = (layer, row, col, cassion_head, linear_CWC)
                    initial_drn.append(drn_input)
            
            if caisson_inactive == True:
                for layer_i in range(idomain.shape[0]):
                    idomain[int(layer_i), int(cassion_row), int(cassion_col)] = 0
        
        origin_drn=[]
        if "drn_0" in self.gwf.package_names:
            drn = self.gwf.get_package('drn_0')
            drn_data = drn.stress_period_data.get_data(kper=0) 
            for i in drn_data[0]:
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
        col = remain_cell % ncol-1
        return layer, row, col
    
    cbc = flopy.utils.CellBudgetFile(f"{gwf.model_ws}/{gwf.name}.cbc")
    drain_leakage = cbc.get_data(text="DRN")
    Drn_out_df = pd.DataFrame({
        'layer': [node_equation(node, nlay, nrow, ncol)[0] for node in drain_leakage[0]["node"]],
        'row': [node_equation(node, nlay, nrow, ncol)[1] for node in drain_leakage[0]["node"]],
        'col': [node_equation(node, nlay, nrow, ncol)[2] for node in drain_leakage[0]["node"]],
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
        self.rw = Cassion_data[0]["well_D"]/2
        self.rs = Cassion_data[0]["rs"]
        self.ab_roughness = Cassion_data[0]["ab_roughness"]
        self.pipe_data_all = Cassion_data[1]
        self.skin_Thick = Cassion_data[0]["skin_Thick"]
        self.pipe_open = Cassion_data[0]["pipe_open"]
        self.aqu_head = aqu_head
        
    def Q_pipe(self):
        def f_calculation(e,pV,D,iq,p):     # Colebrook-White equation
            Re=p*pV*D/u                          # [kg/m^3]*[m/s]*[m]/[kg/m*s] = []
            WRe = p*iq/(math.pi*u)               # [kg/m^3]*[m*2/s]/[kg/m*s] = []
            if pV == 0:
                return 0
            else:
                if Re >= 3000:
                    f1=(1/(-2 * math.log10(e/3.7065 - (5.0272/Re) * math.log10(e/3.827 - 
                       (4.567/Re)*math.log10((e/7.7918)**0.9924 + (5.3326/(208.5 + Re))**0.9345)))))**2
                    f1 = f1*(1-0.0153*WRe**0.3978)
                elif Re <= 2000:
                    f1 = 64/Re*(1+0.04304*WRe**0.6142)
                else:
                    x2=(1/(-2 * math.log10(e/3.7065 - (5.0272/3000) * math.log10(e/3.827 - 
                       (4.567/3000)*math.log10((e/7.7918)**0.9924 + (5.3326/(208.5 + 3000))**0.9345)))))**2
                    x2 = x2*(1-0.0153*WRe**0.3978)
                    x1=64/2000*(1+0.04304*WRe**0.6142)
                    f1=(x2-x1)/(3000-2000)*(Re-2000)+x1
                return f1
        
        def head_loss(IQ, SQ, f2, f1, L, D):
            x = sp.symbols('x')
            f_equation = ((f2 - f1) / L) * x + f1
            Q_equation = (IQ / L) * x + SQ
            Q2_result = sp.integrate(f_equation*Q_equation**2, (x,0,L))

            if friction_loss==True:
                Cf = 8 / (sp.pi ** 2 * g * D ** 5)  
                dhf = (Cf * Q2_result).evalf()
            else:
                dhf=0

            if momemtom_loss==True:
                Cm = 8 / (sp.pi ** 2 * g * D ** 4)
                dhm = (Cm * (IQ**2+2*IQ*SQ)).evalf()
            else:
                dhm=0
            
            return dhf, dhm
                
        for pipe_key, pipe_data in self.pipe_data_all.items():
            cell_count = pipe_data.shape[0]
            for idx, cell in pipe_data.iloc[::-1].iterrows():
                Q_pipe = self.pipe_data_all[pipe_key].loc[idx,"Q_pipe"]
                L = self.pipe_data_all[pipe_key].loc[idx,"length"]
                
                Q_inflow = (self.pipe_data_all[pipe_key].loc[idx,"Q_in"]/L)/86400  #m2/sec
                V_pipe = (self.pipe_data_all[pipe_key].loc[idx,"Q_pipe"] / (math.pi * self.rw ** 2)) / 86400 #([m^3/d]/[m^2])*(d/s)=m/s
                f = f_calculation(self.ab_roughness, V_pipe, self.rw*2, Q_inflow, p)
                
                self.pipe_data_all[pipe_key].loc[idx,"f"] = f2 = f
                
                IQ = self.pipe_data_all[pipe_key].loc[idx,"Q_in"]
                if idx == cell_count-1:
                    f1 = 0
                    SQ = 0
                else:
                    f1 = self.pipe_data_all[pipe_key].loc[idx+1,"f"]
                    SQ = self.pipe_data_all[pipe_key].loc[idx+1,"Q_pipe"]
                    
                dhf, dhm = head_loss(IQ, SQ, f2, f1, L, self.rw*2)
                self.pipe_data_all[pipe_key].loc[idx,"Head_loss_F"] = dhf
                self.pipe_data_all[pipe_key].loc[idx,"Head_loss_M"] = dhm
    
    def Pipe_head_distribute(self):
        for pipe_key, pipe_data in self.pipe_data_all.items():
            V_at_Cassion = (pipe_data.loc[0,"Q_pipe"] /(math.pi * self.rw ** 2))
            Cassion_inflow_head_loss = (V_at_Cassion**2)/(2*g)
            for idx, cell in pipe_data.iterrows():
                row = int(cell["row"])
                col = int(cell["col"])
                layer = int(cell["layer"])
                if idx == 0:
                    if cassion_loss==True:
                        start_head = self.head + Cassion_inflow_head_loss
                    else:
                        start_head = self.head
                else:
                    start_head = self.pipe_data_all[pipe_key].loc[idx-1,"Pipe_head_end"]
                
                dh = self.pipe_data_all[pipe_key].loc[idx,"Head_loss_F"] +self.pipe_data_all[pipe_key].loc[idx,"Head_loss_M"]
                self.pipe_data_all[pipe_key].loc[idx,"Pipe_head_end"] = float(start_head + dh)
                self.pipe_data_all[pipe_key].loc[idx,"Pipe_head_input"] = float(start_head + dh/2)
                
                if self.pipe_data_all[pipe_key].loc[idx,"Pipe_head_input"] > self.aqu_head [layer,row,col]:
                    self.pipe_data_all[pipe_key].loc[idx,"Pipe_head_input"] = float(self.aqu_head [layer,row,col])
                    self.pipe_data_all[pipe_key].loc[idx,"Pipe_head_end"] = float(self.aqu_head [layer,row,col])

    def Pipe_Conductance(self):
        def β_calculate(K):
            K_cs = K*100/86400 #[cm/sec]
            β = 0.5*(K_cs)**(-3/4)*(100/86400)**2  #[day^2/m^2]
            return β
        
        for pipe_key, pipe_data in self.pipe_data_all.items():
            for idx, cell in pipe_data.iterrows():
                K_a = cell["K_gm"]
                β_a = β_calculate(K_a)
                K_s = self.Skin_K
                β_s = β_calculate(K_s)
                rw = self.rw
                rs = self.rs
                re = cell["re"]
                Q_inflow = cell["Q_in"]
                Fe = 1/cell["linear_CWC"]
                L = cell["length"]
                
                if self.skin_Thick==0:
                    rs=self.rw
                    
                if non_linear==True:
                    NS_aq = β_a * (1 / rs - 1 / re)
                    NS_sk = β_s * (1 / rw - 1 / rs)
                else:
                    NS_aq = 0
                    NS_sk = 0
                    
                if slot_loss == True:
                    NS_sl = (1 / (2 * g)) * (1 / (rw * C * self.pipe_open)) ** 2
                else:
                    NS_sl = 0
                    
                common_factor = (1 / (2 * math.pi * L)) ** 2
                Se = Q_inflow*common_factor*(NS_aq + NS_sk + NS_sl)
                CWN=1/(Fe+Se)
                self.pipe_data_all[pipe_key].loc[idx,"S_aq_N"] = Q_inflow*common_factor*NS_aq
                self.pipe_data_all[pipe_key].loc[idx,"S_sk_N"] = Q_inflow*common_factor*NS_sk
                self.pipe_data_all[pipe_key].loc[idx,"S_sl_N"] = Q_inflow*common_factor*NS_sl
                self.pipe_data_all[pipe_key].loc[idx,"CWC"] =CWN

def iterative_data_save(Cassion_input_data, iterative_data, repeat):
    iterative_data[repeat] = {}
    for Cassion_key, Cassion_data in Cassion_input_data.items():
        iterative_data[repeat][Cassion_key] = {}
        for pipe_key, pipe_df in Cassion_data[1].items():
            iterative_data[repeat][Cassion_key][pipe_key] = pipe_df.copy()
            
    return iterative_data
        
        
# ✅ 실행---------------------------------------------------------------------------------------------------------
def run(sim, pipe_excel, sim_name, alpha=alpha, iterative_version=0):

    gwf = sim.get_model(sim_name)
    
    iterative_data = {}
    repeat=0
    old_drn_input_pd = 0
    old_Q_in_max = 0

    root = Radial_well_input(pipe_excel, gwf, sim)
    gwf, origin_drn = root.initial_drn_input()
    success, buff = sim.run_simulation(silent=True)
    if not success:
        print("MODFLOW did not terminate normally.")
    Cassion_input_data, aqu_head, Q_in_max = output_read(root.Cassion_input_data, repeat, gwf)
    iterative_data = iterative_data_save(Cassion_input_data, iterative_data, repeat)
    
    while True:
        old_Q_in_max = Q_in_max
        repeat +=1
        total_pipe_data = pd.DataFrame()
        for cassion_name, Cassion_data in Cassion_input_data.items():
            root = Horizontal_calculate(Cassion_data, aqu_head)
            root.Q_pipe()
            root.Pipe_head_distribute()
            root.Pipe_Conductance()
            Cassion_input_data[cassion_name][1] = root.pipe_data_all
            for pipe_key, pipe_pd in Cassion_input_data[cassion_name][1].items():
                total_pipe_data = pd.concat([total_pipe_data, pipe_pd], ignore_index=True)
        
        iterative_data = iterative_data_save(Cassion_input_data, iterative_data, repeat)
        drn_input_pd = pd.DataFrame(total_pipe_data, columns=["layer", "row", "col", "Pipe_head_input", "CWC"])
        if Damping_start==True and repeat >=2:
            drn_input_pd["Pipe_head_input"] = (1-alpha)*old_drn_input_pd["Pipe_head_input"]+alpha*drn_input_pd["Pipe_head_input"]
        rcw_input = list(drn_input_pd.itertuples(index=False, name=None))
        
        drn_input = origin_drn+rcw_input
        gwf.remove_package("drn_0")
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: drn_input}, save_flows=True)
        sim.write_simulation(silent=True)
        success, buff = sim.run_simulation(silent=True)
        if not success:
            print("MODFLOW did not terminate normally.")
        
        old_drn_input_pd = drn_input_pd
        Cassion_input_data, aqu_head, Q_in_max = output_read(Cassion_input_data, repeat, gwf)
        
        if repeat == 1:
            first_Q_in_max = Q_in_max
            
        if iterative_version==0:
            print(Q_in_max)
            if Q_in_max < convergence_end:
                break
        else:
            Q_max_ratio = Q_in_max/first_Q_in_max
            print(first_Q_in_max, Q_in_max, Q_max_ratio*100)
            if Q_max_ratio < convergence_end_ratio and repeat >=5: 
                break            
            
            
    return iterative_data, Cassion_input_data, sim
        
        
