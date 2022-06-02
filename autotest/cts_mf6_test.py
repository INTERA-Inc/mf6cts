import os
import sys
import platform
import shutil
import string
import time
import numpy as np
import pandas as pd

sys.path.append(".")
import flopy
import pyemu

if "linux" in platform.platform().lower():
    lib_name = os.path.join("..", "bin", "linux", "libmf6.so")
    mt3d_bin = os.path.join("..", "bin", "linux", "mt3dusgs")
    mf6_bin = os.path.join("..", "bin", "linux", "mf6")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    lib_name = os.path.join("..", "bin", "mac", "libmf6.so")
    mt3d_bin = os.path.join("..", "bin", "mac", "mt3dusgs")
    mf6_bin = os.path.join("..", "bin", "mac", "mf6")
else:
    lib_name = os.path.join("..", "bin", "win", "libmf6.dll")
    mt3d_bin = os.path.join("..", "bin", "win", "mt3dusgs.exe")
    mf6_bin = os.path.join("..", "bin", "win", "mf6.exe")

sys.path.append(os.path.join("..", "mf6cts"))

nrowncol = 33
# number of layers
NLAY = 3

# the length of the domain in x and y directions
xlenylen = 10 * nrowncol
# the row and col spacing
delrdelc = xlenylen / nrowncol
# the length of each stress period in days
perlen = [180.0 for _ in range(10)]

# build up the time stepping container for the tdis package
tdis_data = [(p, 1, 1.0) for p in perlen]

# top elev of the model layer 1
top = 10
# total thickness of the domain
z_total = 10

sys_rate_mean = 100.0
np.random.seed(111)
rate_mults = np.random.normal(1.0, 0.25, len(perlen) - 1)
rate_mults[rate_mults < 0] *= -1

gwfname = "gwf"
bud_file = "{}.bud".format(gwfname)
hds_file = "{}.hds".format(gwfname)
gwtname = "gwt"


class Mf6TListBudget(flopy.utils.mflistfile.ListBudget):
    """"""

    def set_budget_key(self):
        self.budgetkey = "MASS BUDGET FOR ENTIRE MODEL"
        return


def setup_five_spotish(plot=False, sim_ws="fivespot", simple_pattern=False, eff1=0.75, eff2=0.95,
                       balance_inout=False, nlay=NLAY, simple_eff=False, ghb_source=False, simple_hk=True,
                       shared_inj=False):
    """create and run a separate flow and tranport sim with
	multiple interior extraction wells and 4 corner injectors.  requires
	the MODFLOW-6 binary to run the model(s)

	"""

    # the number of rows and cols - gotta be odd!

    assert nrowncol % 2 != 0

    # the bottom elev of each layer
    lay_thick = z_total / nlay
    botm = [top - lay_thick]
    for _ in range(1, nlay):
        botm.append(botm[-1] - lay_thick)
    if not simple_hk:
        dxdy = np.zeros(nrowncol) + delrdelc
        gs = pyemu.geostats.GeoStruct(variograms=pyemu.geostats.ExpVario(1.0, delrdelc * 5))
        ss = pyemu.geostats.SpecSim2d(dxdy, dxdy, gs)
        hk = 10 ** ss.draw_arrays(num_reals=nlay, mean_value=1.0)
    else:
        hk = [10 for _ in range(nlay)]

        # if nlay > 1:
        #    hk[1] *= 0.5

    # instaniate...
    sim = flopy.mf6.MFSimulation(sim_name="mfsim", sim_ws=sim_ws, continue_=True, memory_print_option="all")
    tdis = flopy.mf6.ModflowTdis(simulation=sim, nper=len(tdis_data), perioddata=tdis_data)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions="newton")

    # instantiate discretization package

    idomain = np.ones((nlay,nrowncol,nrowncol))
    idomain[:,0,2:5] = 0
    idomain[:,nrowncol-1, 2:5] = 0

    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrowncol, ncol=nrowncol, delr=delrdelc, delc=delrdelc, top=top,
                                  botm=botm[:nlay],idomain=idomain)

    # instantiate node property flow package
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=hk[:nlay], k33overk=True, k33=0.1, icelltype=1, save_specific_discharge=True,
                                  save_flows=True, save_saturation=True)

    # instantiate initial conditions for the flow solution - set starting heads at midpoint of layer 1
    ic = flopy.mf6.ModflowGwfic(gwf, strt=top)

    # instantiate the storage package - stress period 1 is steady state, transient after that...
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1, steady_state={0: True}, transient={1: True}, ss=0.00001, sy=0.01)

    # output control - headsave and budget file names

    oc = flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=bud_file, head_filerecord=hds_file,
                                headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
                                saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
                                printrecord=[("BUDGET", "LAST")], )

    # set up the injection and extraction wells - use the aux var for concentration
    # total in and out well rate in L^3/T

    wel_sp_data = {}
    cts_1, cts_2 = {}, {}
    if simple_pattern:
        for i, rate_mult in enumerate(rate_mults):
            sys_rate = sys_rate_mean
            # if i % 3 == 0:
            #    sys_rate *= 2

            # initialize with the extractors
            # wel_data = [[(nlay - 1, int(nrowncol / 3), int(nrowncol / 3)), -sys_rate, 0.0]]
            # wel_data.append([(nlay - 1, int(nrowncol / 2), int(nrowncol / 2)), -sys_rate, 0.0])
            wel_data = [[(nlay - 1, int(nrowncol / 3), int(2 * nrowncol / 3)), -sys_rate, 0.0]]
            wel_data.append([(nlay - 1, int(2 * nrowncol / 3), int(2 * nrowncol / 3)), -sys_rate, 0.0])


            # # now add the injectors
            if shared_inj:
                wel_data.append([(nlay - 1, 0, int(nrowncol / 2)), sys_rate, 0.0])
                wel_data.append([(nlay - 1, 0, int(nrowncol / 2)), sys_rate, 0.0])

            else:
                if ghb_source:
                    wel_data.append([(nlay - 1, 0, int(nrowncol / 2)), sys_rate, 0.0])
                    wel_data.append([(nlay - 1, nrowncol - 1, int(nrowncol / 2)), sys_rate / 3, 0.0])
                    wel_data.append([(nlay - 1, nrowncol - 1, nrowncol - 1), sys_rate / 3, 0.0])
                    wel_data.append([(nlay - 1, 0, nrowncol - 1), sys_rate / 3, 0.0])
                else:
                    wel_data.append([(nlay - 1, 0, 0), sys_rate, 0.0])
                    wel_data.append([(nlay - 1, nrowncol - 1, 0), sys_rate / 3, 0.0])
                    wel_data.append([(nlay - 1, nrowncol - 1, nrowncol - 1), sys_rate / 3, 0.0])
                    wel_data.append([(nlay - 1, 0, nrowncol - 1), sys_rate / 3, 0.0])

            # if i % 3 == 0:
            #    for iwd,wd in enumerate(wel_data):
            #        wel_data[iwd][1] *= -1

            wel_sp_data[i + 1] = wel_data

            # fill the cts containers - splitting up the ext and inj wells into two systems
            if shared_inj:
                cts_1[i + 2] = [wel_data[0], wel_data[2]]
                cts_2[i + 2] = [wel_data[1], wel_data[3]]
            else:
                cts_1[i + 2] = [wel_data[0], wel_data[2]]
                cts_2[i + 2] = [wel_data[1], wel_data[3],wel_data[4],wel_data[5]]
            # cts_1[i+2] = wel_data
    else:

        for i, rate_mult in enumerate(rate_mults):
            sys_rate = sys_rate_mean * rate_mult
            if not balance_inout and i % 3 == 0:
                sys_rate *= 1

            # initialize with the extractors

            #
            # # now add the injectors
            if balance_inout:
                wel_data = [
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate / 3,
                     0.0]]
                wel_data.append(
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate / 3,
                     0.0])
                wel_data.append(
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate,
                     0.0])
                wel_data.append(
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate / 3,
                     0.0])
                wel_data.append([(0, 0, 0), sys_rate, 0.0])
                wel_data.append([(0, nrowncol - 1, 0), sys_rate / 3., 0.0])
                wel_data.append([(0, nrowncol - 1, nrowncol - 1), sys_rate / 3., 0.0])
                wel_data.append([(0, 0, nrowncol - 1), sys_rate / 3., 0.0])
            else:
                wel_data = [
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate / 2,
                     0.0]]
                wel_data.append(
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate / 2,
                     0.0])
                wel_data.append(
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate / 2,
                     0.0])
                wel_data.append(
                    [(nlay - 1, np.random.randint(0, nrowncol, 1)[0], np.random.randint(0, nrowncol, 1)[0]),
                     -sys_rate / 2,
                     0.0])
                wel_data.append([(0, 0, 0), 2 * sys_rate, 0.0])
                wel_data.append([(0, nrowncol - 1, 0), sys_rate / 3., 0.0])
                wel_data.append([(0, nrowncol - 1, nrowncol - 1), sys_rate / 3., 0.0])
                wel_data.append([(0, 0, nrowncol - 1), sys_rate / 3., 0.0])

            # increment to skip the first stress period
            wel_sp_data[i + 1] = wel_data

            # fill the cts containers - splitting up the ext and inj wells into two systems
            if i % 2 == 0:
                cts_1[i + 2] = [wel_data[0], wel_data[1], wel_data[3], wel_data[4]]
                cts_2[i + 2] = [wel_data[2], wel_data[5], wel_data[6], wel_data[7]]
            else:
                cts_2[i + 2] = [wel_data[0], wel_data[1], wel_data[3], wel_data[4]]
                cts_1[i + 2] = [wel_data[2], wel_data[5], wel_data[6], wel_data[7]]
            # cts_1 = wel_data.copy()
    # instantiate the wel package
    wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data=wel_sp_data, save_flows=True,
                                  auxiliary="concentration", auto_flow_reduce=1.0)

    if ghb_source is not False:
        try:
            gconc = float(ghb_source)
        except:
            gconc = 1.0
        ghb_data = []
        for k in range(nlay):
            # c = 0.0
            # if k == nlay-1:
            #    c = 1.0
            ghb_data.extend([[(k, i, 0), top, 1.0, gconc] for i in range(2, nrowncol - 2)])
            ghb_data.extend([[(k, i, nrowncol - 1), (top + botm[0]) / 2., 1., 0.0] for i in range(2, nrowncol - 2)])
        ghb_perdata = {kper: ghb_data for kper in range(int(len(perlen) / 2))}
        ghb_data = []
        for k in range(nlay):
            ghb_data.extend([[(k, i, 0), top, 5, gconc] for i in range(2, nrowncol - 2)])
            ghb_data.extend([[(k, i, nrowncol - 1), (top + botm[0]) / 2., 5, 0.0] for i in range(2, nrowncol - 2)])
        for kper in range(int(len(perlen) / 2), len(perlen)):
            ghb_perdata[kper] = ghb_data

        ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghb_perdata, auxiliary="concentration", save_flows=True)

    else:
        ghb_data = []
        for k in range(nlay):
            ghb_data.extend([[(k, i, 0), top, 5, 0.0] for i in range(2, nrowncol - 2)])
            ghb_data.extend([[(k, i, nrowncol - 1), (top + botm[0]) / 2., 5, 0.0] for i in range(2, nrowncol - 2)])
        ghb_perdata = {0: ghb_data}

        ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghb_perdata, auxiliary="concentration", save_flows=True)

    # drn_data = [[(0,int(nrowncol/2),int(nrowncol/2)),top,5.0,0.0]]
    # drn = flopy.mf6.ModflowGwfdrn(gwf,stress_period_data=drn_data,auxiliary="concentration",save_flows=True)

    # just a generic solver will do
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab", outer_dvclose=0.0001, inner_dvclose=0.0001,
                               outer_maximum=100,
                               inner_maximum=250)

    # write the input file
    sim.simulation_data.max_columns_of_data = sim.get_model("gwf").dis.nrow.data
    sim.set_all_data_external(check_data=False)
    sim.write_simulation()

    # run the flow model
    # os.chdir(sim_ws)
    shutil.copy2(mf6_bin, os.path.join(sim_ws, os.path.split(mf6_bin)[1]))
    pyemu.os_utils.run("mf6", cwd=sim_ws)
    # os.chdir("..")

    # now make a separate transport simulation and model
    simt_ws = sim_ws + "_t"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    os.makedirs(simt_ws)
    # transport sim
    simt = flopy.mf6.MFSimulation(sim_ws=simt_ws, memory_print_option="all", continue_=True)
    # transport sim temporal discet matches flow solution discret
    tdist = flopy.mf6.ModflowTdis(simulation=simt, nper=len(tdis_data), perioddata=tdis_data)

    # transport model instance
    gwtname = "gwt"
    gwt = flopy.mf6.ModflowGwt(simt, modelname=gwtname, save_flows=True)

    # transport model discret package
    dist = flopy.mf6.ModflowGwtdis(gwt, nlay=nlay, nrow=nrowncol, ncol=nrowncol, delr=delrdelc, delc=delrdelc, top=top,
                                   botm=botm[:nlay],idomain=idomain)

    # copy in the flow solution output files for use in the transport model
    shutil.copy2(os.path.join(sim_ws, hds_file), os.path.join(simt_ws, hds_file))
    shutil.copy2(os.path.join(sim_ws, bud_file), os.path.join(simt_ws, bud_file))
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=[["gwfhead", hds_file], ["gwfbudget", bud_file]],
                                  flow_imbalance_correction=True)

    # initial concen
    strt = 1.0
    if ghb_source:
        strt = 0.0
    ict = flopy.mf6.ModflowGwtic(gwt, strt=strt)

    # remaining transport packages
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.05)
    dsp = flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, alh=1.0, ath1=0.1, ath2=0.1)
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=[["WEL_0", "AUX", "CONCENTRATION"], ["GHB_0", "AUX", "CONCENTRATION"]])

    # transport sim output
    ucn_file = "{}.ucn".format(gwt.name)
    oct = flopy.mf6.ModflowGwtoc(gwt, budget_filerecord="{}.cbc".format(gwt.name),
                                 concentration_filerecord=ucn_file,
                                 concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
                                 saverecord=[("CONCENTRATION", "ALL")],
                                 printrecord=[("BUDGET", "LAST")])

    imst = flopy.mf6.ModflowIms(simt, filename="{}.ims".format(gwt.name), linear_acceleration="bicgstab",
                                inner_dvclose=0.0001, outer_dvclose=0.0001, outer_maximum=100, inner_maximum=100)

    simt.register_ims_package(imst, [gwt.name])

    # write a cts file
    with open(os.path.join(simt_ws, "model.cts"), 'w') as f:
        f.write("begin options\n\nend options\n\n")
        for kper in range(len(perlen) - 1):
            e1 = eff1
            if not simple_eff:
                e1 = min(1, max(0, e1 * np.random.uniform(0, 2.0, 1)[0]))
            f.write("begin period {0} cts 1 efficiency {1:4.3f}\n".format(kper + 2, e1))
            for wd in cts_1[kper + 2]:
                f.write("wel wel_0 {0} {1} {2} {3}\n".format("out" if wd[1] < 0 else "in", wd[0][0] + 1, wd[0][1] + 1,
                                                             wd[0][2] + 1))
            # f.write("{0} {1} {2} {3} DRN_0\n".format(drn_data[0][0][0]+1,drn_data[0][0][1]+1,drn_data[0][0][2]+1,"out"))
            f.write("end period {0} cts 1\n\n".format(kper + 2))

            e2 = eff2
            if not simple_eff:
                e2 = min(1, max(0, e2 * np.random.uniform(0, 2.0, 1)[0]))
            f.write("begin period {0} cts 2 efficiency {1:4.3f}\n".format(kper + 2, e2))
            for wd in cts_2[kper + 2]:
                f.write("wel wel_0 {0} {1} {2} {3}\n".format("out" if wd[1] < 0 else "in", wd[0][0] + 1, wd[0][1] + 1,
                                                             wd[0][2] + 1))
            f.write("end period {0} cts 2\n\n".format(kper + 2))

    # write the transport inputs and run
    simt.write_simulation()
    # os.chdir(simt_ws)
    # os.system("mf6")
    # os.chdir("..")
    shutil.copy2(mf6_bin, os.path.join(simt_ws, os.path.split(mf6_bin)[1]))
    pyemu.os_utils.run("mf6", cwd=simt_ws)


def test_five_spotish_api(prep_pst=False):
    """a basic test of the cts

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    if not prep_pst:
        setup_five_spotish(plot=False, sim_ws=org_sim_ws, eff1=0.5, eff2=0.75, simple_eff=False, nlay=1,
                           ghb_source=True,
                           simple_pattern=True, simple_hk=False)
    else:
        setup_five_spotish(plot=False, sim_ws=org_sim_ws, eff1=0.7,eff2=0.7, simple_eff=True, nlay=1,
                           ghb_source=100,
                           simple_pattern=True, simple_hk=True)

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = org_sim_ws + "_t"
    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0, axis=1),
                   :].mass.sum()
    out_node_mass = node_df.loc[node_df.apply(lambda x: x.cum_vol < 0, axis=1),
                    :].mass.sum()
    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    if in_node_mass != 0.0:
        assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    eff = sys_df.mass_treated / sys_df.mass_removed
    d = (eff - sys_df.requested_efficiency).apply(np.abs)
    print(d)
    assert d.max() < 1.0e-10

    d = (sys_df.mass_removed - (sys_df.mass_treated + sys_df.mass_injected)).apply(np.abs)
    print(d)

    wel_df = pd.read_csv(os.path.join(sim_ws,"gwf.wel_stress_period_data_2.txt"),header=None,delim_whitespace=True,names=["l","r","c","flux","concen"])
    print(wel_df)

    node_df = pd.read_csv(os.path.join(sim_ws, "gwf_cts_flow_node_summary.csv"))
    cts_vals = node_df.loc[node_df.stress_period==2,"requested_rate"]
    cts_vals.values.sort()
    wel_vals = wel_df.flux.values
    wel_vals.sort()
    d = np.abs(wel_vals - cts_vals).sum()
    print(d)
    assert d < 1.0e-10

def plot(simt_ws, sim_ws):
    """ a plotting function

    """
    wel_data = {}
    if "maw" not in simt_ws.lower():
        with open(os.path.join(sim_ws, "gwf.wel"), 'r') as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                if line.lower().strip().startswith("begin period"):
                    wd = []
                    period = int(line.strip().split()[-1])
                    while True:
                        line2 = f.readline()
                        if line2 == "":
                            raise Exception()
                        if line2.lower().strip().startswith("end period"):
                            break
                        raw = line2.strip().split()
                        l, r, c = int(raw[0]) - 1, int(raw[1]) - 1, int(raw[2]) - 1
                        flx = float(raw[3])
                        wd.append([l, r, c, flx])
                    wel_data[period] = wd
    else:
        with open(os.path.join(sim_ws, "gwf.maw"), 'r') as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                if line.lower().strip().startswith("begin connectiondata"):
                    wd = []
                    while True:
                        line2 = f.readline()
                        if line2 == "":
                            raise Exception()
                        if line2.lower().strip().startswith("end connectiondata"):
                            break
                        raw = line2.strip().split()
                        l, r, c = int(raw[2]) - 1, int(raw[3]) - 1, int(raw[4]) - 1
                        wd.append([l, r, c, 1])
                    break
        for i in range(1000):
            wel_data[i] = wd

    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    # some cheap plotting...
    hds = flopy.utils.HeadFile(os.path.join(simt_ws, "gwf.hds"), precision="double")
    ucn = flopy.utils.HeadFile(os.path.join(simt_ws, "gwt.ucn"), precision="double", text="concentration")

    with PdfPages(os.path.join(simt_ws, "results.pdf")) as pdf:
        hds_times = hds.get_times()
        ucn_times = ucn.get_times()
        hmin, hmax = 1.0e+20, -1.0e+20
        cmin, cmax = 1.0e+20, -1.0e+20
        for htime, ctime in zip(hds_times, ucn_times):
            hmin = min(hmin, hds.get_data(totim=htime).min())
            hmax = max(hmax, hds.get_data(totim=htime).max())
            cmin = min(cmin, ucn.get_data(totim=ctime).min())
            cmax = max(cmax, ucn.get_data(totim=ctime).max())

        # hlevels = np.linspace(hmin, hmax, 10)
        hlevels = 4
        for itime, (htime, ctime) in enumerate(zip(hds_times, ucn_times)):
            print(htime, ctime)
            wd = wel_data.get(itime + 1, [])

            fig, axes = plt.subplots(hds.get_data().shape[0], 1, figsize=(7, 7 * hds.get_data().shape[0]))
            axes = np.atleast_1d(axes)
            hdata = hds.get_data(totim=htime)
            cdata = ucn.get_data(totim=ctime)
            for k in range(hdata.shape[0]):
                wel_mask = np.zeros((nrowncol, nrowncol))
                for item in wd:
                    if item[0] != k:
                        continue
                    wel_mask[item[1], item[2]] = item[3]
                wel_mask[wel_mask == 0] = np.NaN
                cb = axes[k].imshow(cdata[k, :, :], vmax=cmax, vmin=cmin, cmap="cool")
                cb = plt.colorbar(cb, ax=axes[k])
                cb.set_label("concentration")
                axes[k].set_title("concentration at {0} layer {1}".format(ctime, k + 1))
                axes[k].imshow(wel_mask, cmap="jet")
                cs = axes[k].contour(hdata[k, :, :], levels=hlevels, colors="k")
                axes[k].clabel(cs, cs.levels)

            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)


def convert_5spot_to_maw(org_sim_ws="fivespot", nlay=NLAY, eff1=0.0, eff2=0.0, simple_eff=True):
    """conver from WEL to MAW

    """
    # name of the flow solution files
    gwfname = "gwf"

    sim_ws = org_sim_ws + "_maw"
    org_simt_ws = org_sim_ws + "_t"

    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)

    # instaniate...
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws)
    gwf = sim.get_model(gwfname)

    # wdata = gwf.get_package("WEL_0").stress_period_data
    gwf.remove_package("WEL_0")

    lay_thick = z_total / nlay
    botm = [top - lay_thick]
    for _ in range(1, nlay):
        botm.append(botm[-1] - lay_thick)
    hk = [10 for _ in range(nlay)]
    if nlay > 1:
        hk[1] = 0.1

    # build up maw inputs
    maw_pd = []
    maw_cd = []
    inj_is = [0, nrowncol - 1, nrowncol - 1]
    inj_js = [0, nrowncol - 1, nrowncol - 1, 0]
    inj_bnames = ["cts1", "cts2", "cts2", "cts2"]
    # inj_bnames = ["cts1"]#, "cts1", "cts1", "cts1"]
    imaw = 0
    cts_lines = {"cts1": [], "cts2": []}
    inj_wellno = []
    mawgwfname = "MAW_0"
    mawbudfile = gwfname + ".maw.bud"

    for i, j, bname in zip(inj_is, inj_js, inj_bnames):
        # <wellno> <radius> <bottom> <strt> <condeqn> <ngwfnodes> aux-concen
        maw_pd.append([imaw, delrdelc / 4, botm[-1], top, "specified", nlay, 0.0, bname])
        for k in range(nlay):
            # wellno icon lay row col scrn_top scrn_bot hk_skin radius_skin
            maw_cd.append([imaw, k, (k, i, j), -1e+30, -1e+30, hk[k], 1.0e+30])
        inj_wellno.append(imaw)
        cts_lines[bname].append("maw {0} in {1}\n".format(mawgwfname, imaw + 1))
        imaw += 1
    ext_is = [int(nrowncol / 4), int(7 * nrowncol / 8)]
    ext_js = [int(nrowncol / 4), int(nrowncol / 2)]
    ext_bnames = ["cts1", "cts2"]
    ext_wellno = []
    for i, j, bname in zip(ext_is, ext_js, ext_bnames):
        # <wellno> <radius> <bottom> <strt> <condeqn> <ngwfnodes> aux-concen
        maw_pd.append([imaw, delrdelc / 4, botm[-1], top, "specified", nlay, 1.0e+30, bname])
        for k in range(nlay):
            # wellno icon lay row col scrn_top scrn_bot hk_skin radius_skin
            maw_cd.append([imaw, k, (k, i, j), -1e+30, -1e+30, hk[k], 1.0e+30])
        ext_wellno.append(imaw)
        cts_lines[bname].append("maw {0} out {1}\n".format(mawgwfname, imaw + 1))
        imaw += 1

    maw_perdata = {}
    num_inj = len(inj_wellno)
    cts_kper_dict = {}
    for iper, rmult in enumerate(rate_mults):
        perdata = []
        ctsdata = []
        for wellno in inj_wellno:
            if iper == 0:
                perdata.append([wellno, "status", "inactive"])
            else:
                perdata.append([wellno, "rate", sys_rate_mean * rmult / num_inj])
                perdata.append([wellno, "status", "active"])
                ctsdata.append("maw {0} in {1}\n".format(mawgwfname, wellno))
        for wellno in ext_wellno:
            if iper == 0:
                perdata.append([wellno, "status", "inactive"])
            else:
                perdata.append([wellno, "rate", -2. * sys_rate_mean * rmult])
                perdata.append([wellno, "status", "active"])
                ctsdata.append("maw {0} out {1}\n".format(mawgwfname, wellno))
        maw_perdata[iper] = perdata
        cts_kper_dict[iper] = ctsdata

    flopy.mf6.ModflowGwfmaw(gwf, auxiliary="concentration", save_flows=True,
                            packagedata=maw_pd, connectiondata=maw_cd,
                            perioddata=maw_perdata, boundnames=True, budget_filerecord=mawbudfile, pname=mawgwfname)

    # return
    # write the input file
    sim.write_simulation()

    # run the flow model
    # os.chdir(sim_ws)
    # os.system("mf6")
    # os.chdir("..")
    shutil.copy2(mf6_bin, os.path.join(sim_ws, os.path.split(mf6_bin)[1]))
    pyemu.os_utils.run("mf6", cwd=sim_ws)

    simt_ws = sim_ws + "_t"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)

    shutil.copytree(org_simt_ws, simt_ws)

    shutil.copy2(os.path.join(sim_ws, hds_file), os.path.join(simt_ws, hds_file))
    shutil.copy2(os.path.join(sim_ws, bud_file), os.path.join(simt_ws, bud_file))
    shutil.copy2(os.path.join(sim_ws, mawbudfile), os.path.join(simt_ws, mawbudfile))

    sim = flopy.mf6.MFSimulation.load(sim_ws=simt_ws)
    gwt = sim.get_model(gwtname)
    gwt.remove_package("SSM")
    flopy.mf6.ModflowGwtssm(gwt, sources=[["GHB_0", "AUX", "CONCENTRATION"]])
    gwt.remove_package("FMI")
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=[["gwfhead", hds_file], ["gwfbudget", bud_file],
                                                    [mawgwfname, mawbudfile]],
                                  flow_imbalance_correction=True)

    mwt_pakdata = [[wn, 0.0] for wn in ext_wellno]
    mwt_pakdata.extend([[wn, 0.0] for wn in inj_wellno])

    mwt_perdata = {}
    for ikper in range(len(perlen)):
        items = []
        for wn in inj_wellno:
            if ikper == 0:
                items.append([wn, "status", "inactive"])
            else:
                items.append([wn, "status", "active"])
                items.append([wn, "concentration", 0.0])
                items.append([wn, "rate", 0.0])

        mwt_perdata[ikper] = items
    flopy.mf6.ModflowGwtmwt(gwt, packagedata=mwt_pakdata, mwtperioddata=mwt_perdata, flow_package_name=mawgwfname)

    cts_names = list(cts_lines.keys())
    cts_names.sort()
    with open(os.path.join(simt_ws, "model.cts"), 'w') as f:
        f.write("begin options\n\nend options\n\n")
        for ikper in range(1, len(perlen)):
            e1 = eff1
            e2 = eff2
            if not simple_eff:
                e1 *= np.random.uniform(0.1, 2.0)
                e1 = max(e1, 0)
                e1 = min(e1, 1)
                e2 *= np.random.uniform(0.1, 2.0)
                e2 = max(e2, 0)
                e2 = min(e2, 1)

            for cts_name, eff in zip(cts_names, [e1, e2]):
                f.write("begin period {0} cts {1} efficiency {2}\n".format(ikper + 1, cts_name[-1], eff))
                for item in cts_lines[cts_name]:
                    f.write("  " + item)
                f.write("end period {0}\n\n".format(ikper + 1))

    sim.write_simulation()
    # os.chdir(simt_ws)
    # os.system("mf6")
    # os.chdir("..")
    shutil.copy2(mf6_bin, os.path.join(simt_ws, os.path.split(mf6_bin)[1]))
    pyemu.os_utils.run("mf6", cwd=simt_ws)


def test_five_spotish_api_maw():
    """Test the MAW formulation


    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, nlay=1)
    convert_5spot_to_maw("fivespot", nlay=1)
    org_sim_ws = org_sim_ws + "_maw"
    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = org_sim_ws + "_t"
    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    shutil.copy2(os.path.join(sim_ws, "gwf.maw.bud"), os.path.join(simt_ws, "gwf.maw.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    # assert np.abs(cum.loc[cum.index[-1], "maw"]) > np.abs(api_cum.loc[api_cum.index[-1], "maw"])
    abs_fr_diff = np.abs(api_cum.loc[api_cum.index[-1], "maw"]) / np.abs(cum.loc[cum.index[-1], "maw"])
    assert abs_fr_diff < 0.025, abs_fr_diff

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01, abs_frac_diff

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01


def test_five_spotish_simple_api1():
    """More advanced WEL usecase testing

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.dll"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, simple_pattern=True, eff1=0.0, eff2=0.0, nlay=1)

    t = org_sim_ws + "simple"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_sim_ws, t)
    base_sim_ws = org_sim_ws
    org_sim_ws = t

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = base_sim_ws + "_t"
    t = base_sim_ws + "simple_t"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_simt_ws, t)
    org_simt_ws = t
    t = None

    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()

    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    assert np.abs(cum.loc[cum.index[-1], "wel"]) > np.abs(api_cum.loc[api_cum.index[-1], "wel"])
    abs_frac_diff = np.abs(api_cum.loc[api_cum.index[-1], "wel"] / cum.loc[cum.index[-1], "wel"])
    assert abs_frac_diff < 0.01

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01, abs_frac_diff

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    # now reset with some


def test_five_spotish_simple_api2():
    """Still more WEL use case testing, including "reusing" cts info across stress periods


    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, simple_pattern=True, eff1=0.5, eff2=0.5, simple_eff=True, nlay=1)

    t = org_sim_ws + "simple"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_sim_ws, t)
    base_sim_ws = org_sim_ws
    org_sim_ws = t

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = base_sim_ws + "_t"
    t = base_sim_ws + "simple_t"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_simt_ws, t)
    org_simt_ws = t
    t = None

    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    assert np.abs(cum.loc[cum.index[-1], "wel"]) > np.abs(api_cum.loc[api_cum.index[-1], "wel"])
    abs_frac_diff = np.abs(api_cum.loc[api_cum.index[-1], "wel"] / cum.loc[cum.index[-1], "wel"])
    assert abs_frac_diff < 0.01

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert np.isclose(abs_frac_diff, 1.0)

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert np.isclose(abs_frac_diff, 1.0)

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    # rewrite the cts file to test reuse
    temp_ws = simt_ws + "_temp"
    if os.path.exists(temp_ws):
        shutil.rmtree(temp_ws)
    shutil.copytree(simt_ws, temp_ws)
    cts_file = os.path.join(simt_ws, "model.cts")
    lines = open(cts_file, 'r').readlines()
    with open(cts_file, 'w') as f:
        for line in lines:
            f.write(line)
            if "end period 2 cts 2" in line.lower():
                break

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None

    org_node_df = pd.read_csv(os.path.join(temp_ws, "gwt_cts_node_summary.csv"))
    new_node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))

    d = (org_node_df.iloc[:, 7:] - new_node_df.iloc[:, 7:]).apply(np.abs)

    print(d.max())
    assert d.max().max() < 1.0e-6, d.max().max()


def test_five_spotish_simple_api_mk2k_compare():
    """A version of the simple mf2k-mst test

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, simple_pattern=True, eff1=0.5, eff2=0.5, simple_eff=True, nlay=1)

    t = org_sim_ws + "simple"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_sim_ws, t)
    base_sim_ws = org_sim_ws
    org_sim_ws = t

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = base_sim_ws + "_t"
    t = base_sim_ws + "simple_t"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_simt_ws, t)
    org_simt_ws = t
    t = None

    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    assert np.abs(cum.loc[cum.index[-1], "wel"]) > np.abs(api_cum.loc[api_cum.index[-1], "wel"])
    abs_frac_diff = np.abs(api_cum.loc[api_cum.index[-1], "wel"] / cum.loc[cum.index[-1], "wel"])
    assert abs_frac_diff < 0.01

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert np.isclose(abs_frac_diff, 1.0)

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert np.isclose(abs_frac_diff, 1.0)

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01


def invest():
    conc = flopy.utils.HeadFile(os.path.join("fivespotsimple_t_api", "gwt.ucn"), text="concentration")
    for time in conc.get_times():
        arr = conc.get_data(totim=time)
        print(time, arr[2, 17, 17], arr[0, 0, 0])


def test_five_spotish_api_maw_configfile():
    """Test MAW with the config file usage

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, nlay=1)
    convert_5spot_to_maw("fivespot", nlay=1)
    org_sim_ws = org_sim_ws + "_maw"
    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = org_sim_ws + "_t"
    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    # first copy the mf6cts.py file to this dir
    src_file = os.path.join("..", "mf6cts", "mf6cts.py")
    dest_file = "mf6cts.py"
    if os.path.exists(dest_file):
        os.remove(dest_file)
    shutil.copy2(src_file, dest_file)

    # now write a config file
    config_file = "config_file.py"
    if os.path.exists(config_file):
        os.remove(config_file)
    with open(config_file, 'w') as f:
        f.write("cts_filename='{0}'\n".format("model.cts"))
        f.write("lib_name='{0}'\n".format(os.path.split(lib_name)[-1]))
        f.write("transport_dir='{0}'\n".format(simt_ws))
        f.write("flow_dir='{0}'\n".format(sim_ws))
        f.write("is_structured=True\n")
        f.write("flow_output_files=['gwf.hds','gwf.bud','gwf.maw.bud']\n")

    for fname in ['gwf.hds', 'gwf.bud', 'gwf.maw.bud']:
        if os.path.exists(os.path.join(simt_ws, fname)):
            os.remove(os.path.join(simt_ws, fname))

    os.system("python mf6cts.py config_file.py")
    os.remove(dest_file)

    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    # assert np.abs(cum.loc[cum.index[-1], "maw"]) > np.abs(api_cum.loc[api_cum.index[-1], "maw"])
    abs_fr_diff = np.abs(api_cum.loc[api_cum.index[-1], "maw"]) / np.abs(cum.loc[cum.index[-1], "maw"])
    assert abs_fr_diff < 0.025, abs_fr_diff

    node_df = pd.read_csv(os.path.join(sim_ws, "gwf_cts_flow_node_summary.csv"))
    node_df = node_df.loc[node_df.requested_rate > 0, :]
    d = (node_df.requested_rate - node_df.actual_rate).apply(np.abs).sum()
    print(d)
    assert d > 1000

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01, abs_frac_diff

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01


def test_five_spotish_api_maw_configfile_unbalanced():
    """Test MAW with the config file usage

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, nlay=1)
    convert_5spot_to_maw("fivespot", nlay=1)
    org_sim_ws = org_sim_ws + "_maw"
    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = org_sim_ws + "_t"
    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    # first copy the mf6cts.py file to this dir
    src_file = os.path.join("..", "mf6cts", "mf6cts.py")
    dest_file = "mf6cts.py"
    if os.path.exists(dest_file):
        os.remove(dest_file)
    shutil.copy2(src_file, dest_file)

    # now write a config file
    config_file = "config_file.py"
    if os.path.exists(config_file):
        os.remove(config_file)
    with open(config_file, 'w') as f:
        f.write("cts_filename='{0}'\n".format("model.cts"))
        f.write("lib_name='{0}'\n".format(os.path.split(lib_name)[-1]))
        f.write("transport_dir='{0}'\n".format(simt_ws))
        f.write("flow_dir='{0}'\n".format(sim_ws))
        f.write("is_structured=True\n")
        f.write("flow_output_files=['gwf.hds','gwf.bud','gwf.maw.bud']\n")
        f.write("balance_flows=False\n")

    for fname in ['gwf.hds', 'gwf.bud', 'gwf.maw.bud']:
        if os.path.exists(os.path.join(simt_ws, fname)):
            os.remove(os.path.join(simt_ws, fname))

    os.system("python mf6cts.py config_file.py")
    os.remove(dest_file)

    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    node_df = pd.read_csv(os.path.join(sim_ws, "gwf_cts_flow_node_summary.csv"))
    node_df = node_df.loc[node_df.requested_rate > 0,:]
    d = (node_df.requested_rate - node_df.actual_rate).apply(np.abs).sum()
    print(d)
    assert d < 1.0e-10

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    # assert np.abs(cum.loc[cum.index[-1], "maw"]) > np.abs(api_cum.loc[api_cum.index[-1], "maw"])
    abs_fr_diff = np.abs(api_cum.loc[api_cum.index[-1], "maw"]) / np.abs(cum.loc[cum.index[-1], "maw"])

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)


    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)

def test_five_spotish_complex_api_mk2k_compare():
    """more complex mf2k-mst test

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, simple_pattern=False,
                       eff1=0.95, eff2=0.75, nlay=3, balance_inout=True)

    t = org_sim_ws + "simple"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_sim_ws, t)
    base_sim_ws = org_sim_ws
    org_sim_ws = t

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = base_sim_ws + "_t"
    t = base_sim_ws + "simple_t"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_simt_ws, t)
    org_simt_ws = t
    t = None

    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    assert np.isclose(np.abs(api_cum.loc[api_cum.index[-1], "wel"]), 0.0)
    # abs_frac_diff = np.abs(api_cum.loc[api_cum.index[-1], "wel"]/cum.loc[cum.index[-1], "wel"])
    # assert abs_frac_diff < 0.01,abs_frac_diff

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)

    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01


def test_five_spotish_api_maw_complex():
    """More complex MAW use case testing including checking the maw-specific node summary file

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, nlay=3, ghb_source=True)
    convert_5spot_to_maw("fivespot", nlay=3, eff1=0.5, eff2=0.95, simple_eff=False)
    org_sim_ws = org_sim_ws + "_maw"
    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = org_sim_ws + "_t"
    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    shutil.copy2(os.path.join(sim_ws, "gwf.maw.bud"), os.path.join(simt_ws, "gwf.maw.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    # assert np.abs(cum.loc[cum.index[-1], "maw"]) > np.abs(api_cum.loc[api_cum.index[-1], "maw"])
    abs_fr_diff = np.abs(api_cum.loc[api_cum.index[-1], "maw"]) / np.abs(cum.loc[cum.index[-1], "maw"])
    assert abs_fr_diff < 0.02, abs_fr_diff

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.02

    # abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_OUT.max()) / out_node_mass)
    # print(abs_frac_diff)
    # assert abs_frac_diff < 0.1

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.02

    eff = sys_df.mass_treated / sys_df.mass_removed
    d = (eff - sys_df.requested_efficiency).apply(np.abs)
    print(d)
    assert d.max() < 1.0e-10

    d = (sys_df.mass_removed - (sys_df.mass_treated + sys_df.mass_injected)).apply(np.abs)
    print(d)
    assert d.max() < 1.0e-10

    maw_df = pd.read_csv(os.path.join(simt_ws, "gwt_maw_node_summary.csv"))

    for wellno in node_df.loc[:, "index"]:
        wdf = maw_df.loc[maw_df.maw_wellno == wellno]
        wdf = wdf.loc[wdf.stress_period == wdf.stress_period.max()]
        wdf = wdf.sum()
        ndf = node_df.loc[node_df.loc[:, "index"] == wellno, :]
        ndf = ndf.loc[ndf.stress_period == ndf.stress_period.max(), :]
        d = np.abs(wdf.cum_mass - ndf.cum_mass).values[0]
        print(wellno, d)
        assert d < 1.0e-10, d
        # print(wdf)
        # print(ndf)


def plot_domain(sim_ws):
    """ a useless function

    """
    import matplotlib.pyplot as plt
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws)
    gwf = sim.get_model("gwf")
    hk = np.log10(gwf.get_package("npf").k.array)
    mn, mx = hk.min(), hk.max()
    print(hk.shape)
    fig, axes = plt.subplots(1, 3, figsize=(7, 2))
    ghb_arr = np.zeros((hk.shape[1], hk.shape[2]))
    ghb_arr[1:-1, [0, -1]] = 1
    ghb_arr = np.ma.masked_where(ghb_arr == 0, ghb_arr)
    hds = flopy.utils.HeadFile(os.path.join(sim_ws, "gwf.hds"), precision="double")
    harr = hds.get_data((0, 0))
    for ax, hkk, aval, ha in zip(axes, hk, ["A) layer 1", "B) layer 2", "C) layer 3"], harr):
        cb = ax.imshow(hkk, vmin=mn, vmax=mx)
        ax.imshow(ghb_arr, cmap="cool_r")
        ax.set_title(aval, loc="left")
        ax.set_xlabel("row")
        ax.set_ylabel("column")
        ct = ax.contour(ha, levels=4, colors="k", linewidths=0.5)
        ax.clabel(ct)

    plt.colorbar(cb, ax=ax, label="$log_{10}\\frac{m}{d}$")
    plt.tight_layout()
    plt.savefig("domain.pdf")


def plot_results_pub(simt_ws, sim_ws):
    """more plotting

    """
    wel_data = {}
    if "maw" not in simt_ws.lower():
        with open(os.path.join(sim_ws, "gwf.wel"), 'r') as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                if line.lower().strip().startswith("begin period"):
                    wd = []
                    period = int(line.strip().split()[-1])
                    while True:
                        line2 = f.readline()
                        if line2 == "":
                            raise Exception()
                        if line2.lower().strip().startswith("end period"):
                            break
                        raw = line2.strip().split()
                        l, r, c = int(raw[0]) - 1, int(raw[1]) - 1, int(raw[2]) - 1
                        flx = float(raw[3])
                        wd.append([l, r, c, flx])
                    wel_data[period] = wd
    else:
        with open(os.path.join(sim_ws, "gwf.maw"), 'r') as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                if line.lower().strip().startswith("begin connectiondata"):
                    wd = []
                    while True:
                        line2 = f.readline()
                        if line2 == "":
                            raise Exception()
                        if line2.lower().strip().startswith("end connectiondata"):
                            break
                        raw = line2.strip().split()
                        l, r, c = int(raw[2]) - 1, int(raw[3]) - 1, int(raw[4]) - 1
                        wd.append([l, r, c, 1])
                    break
        for i in range(1000):
            wel_data[i] = wd

    import string
    import matplotlib.pyplot as plt
    # some cheap plotting...
    hds = flopy.utils.HeadFile(os.path.join(simt_ws, "gwf.hds"), precision="double")
    ucn = flopy.utils.HeadFile(os.path.join(simt_ws, "gwt.ucn"), precision="double", text="concentration")
    hds_times = hds.get_times()[::4]
    ucn_times = ucn.get_times()[::4]
    hmin, hmax = 1.0e+20, -1.0e+20
    cmin, cmax = 1.0e+20, -1.0e+20
    for htime, ctime in zip(hds_times, ucn_times):
        hmin = min(hmin, hds.get_data(totim=htime).min())
        hmax = max(hmax, hds.get_data(totim=htime).max())
        cmin = min(cmin, ucn.get_data(totim=ctime).min())
        cmax = max(cmax, ucn.get_data(totim=ctime).max())
    fig, axes = plt.subplots(len(hds_times), hds.get_data().shape[0], figsize=(7.5, 6))
    axes = np.atleast_2d(axes).T
    hlevels = 4
    acount = 0

    abet = string.ascii_uppercase
    for itime, (htime, ctime) in enumerate(zip(hds_times, ucn_times)):
        print(htime, ctime)
        wd = wel_data.get(itime + 1, [])
        hdata = hds.get_data(totim=htime)
        cdata = ucn.get_data(totim=ctime)
        for k in range(hdata.shape[0]):
            wel_mask = np.zeros((nrowncol, nrowncol)) - 1
            for item in wd:
                if item[0] != k:
                    continue
                if item[3] > 0:
                    wel_mask[item[1], item[2]] = 1
                else:
                    wel_mask[item[1], item[2]] = 0
            wel_mask = np.ma.masked_where(wel_mask == -1, wel_mask)
            cb = axes[itime, k].imshow(cdata[k, :, :], vmax=cmax, vmin=cmin, cmap="cividis", alpha=0.85)
            cb = plt.colorbar(cb, ax=axes[itime, k])
            cb.set_label("concentration")
            axes[itime, k].set_title("{2}) time:{0} days, layer:{1}".format(ctime, k + 1, abet[acount]), loc="left")
            axes[itime, k].imshow(wel_mask, cmap="cool", vmin=0, vmax=1)
            cs = axes[itime, k].contour(hdata[k, :, :], levels=hlevels, colors="k", linewidths=0.25)
            axes[itime, k].clabel(cs, cs.levels)
            axes[itime, k].set_xlabel("column")
            axes[itime, k].set_ylabel("row")

            acount += 1

    plt.tight_layout()
    plt.savefig(os.path.join(simt_ws, "results_pub.pdf"))


def run():
    """function to run balanced and unbalance flows within the pest framework

    """
    class Mf6TListBudget(flopy.utils.mflistfile.ListBudget):
        """"""

        def set_budget_key(self):
            self.budgetkey = "MASS BUDGET FOR ENTIRE MODEL"
            return

    simt_ws = "fivespot_t_api"
    sim_ws = "fivespot_api"

    pyemu.os_utils.run("python mf6cts.py config_file.py")
    hds = flopy.utils.HeadFile(os.path.join(simt_ws, "gwf.hds"), precision="double")
    ucn = flopy.utils.HeadFile(os.path.join(simt_ws, "gwt.ucn"), precision="double", text="concentration")
    arr = hds.get_data()[0]
    np.savetxt(os.path.join("apiheads.dat"), arr)

    arr = ucn.get_data()[0]
    np.savetxt(os.path.join("apiconcen.dat"), arr)

    lst = Mf6TListBudget(os.path.join(simt_ws, "gwt.lst"))
    inc, cum = lst.get_dataframes(start_datetime=None)
    inc.index.name = "totim"
    cum.index.name = "totim"
    inc.columns = inc.columns.map(lambda x: x.lower().replace("_", "-"))
    inc.to_csv("apigwtinc.csv")
    cum.columns = cum.columns.map(lambda x: x.lower().replace("_", "-"))
    cum.to_csv("apigwtcum.csv")

    lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, "gwf.lst"))
    inc, cum = lst.get_dataframes(start_datetime=None)
    inc.index.name = "totim"
    cum.index.name = "totim"
    inc.columns = inc.columns.map(lambda x: x.lower().replace("_", "-"))
    inc.to_csv("apigwfinc.csv")
    cum.columns = cum.columns.map(lambda x: x.lower().replace("_", "-"))
    cum.to_csv("apigwfcum.csv")

    for ws in [sim_ws, simt_ws]:
        csv_files = [f for f in os.listdir(ws) if f.lower().endswith(".csv")]
        for csv_file in csv_files:
            shutil.copy2(os.path.join(ws, csv_file), csv_file)

    # pyemu.os_utils.run("mf6",cwd=sim_ws)
    pyemu.os_utils.run("python mf6cts.py config_file_unbalanced.py")
    # flow_output_files = ['gwf.hds', 'gwf.bud']
    # for f in flow_output_files:
    #     shutil.copy2(os.path.join(sim_ws,f),os.path.join(simt_ws,f))
    # pyemu.os_utils.run("mf6", cwd=simt_ws)
    hds = flopy.utils.HeadFile(os.path.join(simt_ws, "gwf.hds"), precision="double")
    ucn = flopy.utils.HeadFile(os.path.join(simt_ws, "gwt.ucn"), precision="double", text="concentration")
    arr = hds.get_data()[0]
    np.savetxt(os.path.join("heads.dat"), arr)

    arr = ucn.get_data()[0]
    np.savetxt(os.path.join("concen.dat"), arr)

    lst = Mf6TListBudget(os.path.join(simt_ws, "gwt.lst"))
    inc, cum = lst.get_dataframes(start_datetime=None)
    inc.index.name = "totim"
    cum.index.name = "totim"
    inc.columns = inc.columns.map(lambda x: x.lower().replace("_", "-"))
    inc.to_csv("gwtinc.csv")
    cum.columns = cum.columns.map(lambda x: x.lower().replace("_", "-"))
    cum.to_csv("gwtcum.csv")

    lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, "gwf.lst"))
    inc, cum = lst.get_dataframes(start_datetime=None)
    inc.index.name = "totim"
    cum.index.name = "totim"
    inc.columns = inc.columns.map(lambda x: x.lower().replace("_", "-"))
    inc.to_csv("gwfinc.csv")
    cum.columns = cum.columns.map(lambda x: x.lower().replace("_", "-"))
    cum.to_csv("gwfcum.csv")


def invest_run(t_d):
    """test the run() function
    """

    b_d = os.getcwd()
    os.chdir(t_d)
    run()
    os.chdir(b_d)


def setup_pst():
    """setup the pest interface for the prior monte carlo
    """
    test_five_spotish_api(prep_pst=True)
    org_sim_ws = "fivespot_api"
    org_simt_ws = "fivespot_t_api"

    t_d = "temp"
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    os.mkdir(t_d)
    sim_ws = os.path.join(t_d, "fivespot_api")
    simt_ws = os.path.join(t_d, "fivespot_t_api")
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copytree(org_simt_ws, simt_ws)

    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws)

    sim.set_all_data_external(check_data=False)
    sim.simulation_data.max_columns_of_data = sim.get_model("gwf").dis.nrow.data
    sim.write_simulation()

    sim = flopy.mf6.MFSimulation.load(sim_ws=simt_ws)
    sim.set_all_data_external(check_data=False)
    sim.write_simulation()

    config_file = os.path.join(t_d, "config_file.py")
    if os.path.exists(config_file):
        os.remove(config_file)
    with open(config_file, 'w') as f:
        f.write("cts_filename='{0}'\n".format("model.cts"))
        f.write("lib_name='{0}'\n".format(os.path.split(lib_name)[-1]))
        f.write("transport_dir='{0}'\n".format(os.path.split(simt_ws)[-1]))
        f.write("flow_dir='{0}'\n".format(os.path.split(sim_ws)[-1]))
        f.write("is_structured=True\n")
        f.write("flow_output_files=['gwf.hds','gwf.bud']\n")
    config_file = os.path.join(t_d, "config_file_unbalanced.py")
    if os.path.exists(config_file):
        os.remove(config_file)
    with open(config_file, 'w') as f:
        f.write("cts_filename='{0}'\n".format("model.cts"))
        f.write("lib_name='{0}'\n".format(os.path.split(lib_name)[-1]))
        f.write("transport_dir='{0}'\n".format(os.path.split(simt_ws)[-1]))
        f.write("flow_dir='{0}'\n".format(os.path.split(sim_ws)[-1]))
        f.write("is_structured=True\n")
        f.write("flow_output_files=['gwf.hds','gwf.bud']\n")
        f.write("balance_flows=False\n ")
    shutil.copy2(os.path.join("..", "mf6cts", "mf6cts.py"), os.path.join(t_d, "mf6cts.py"))
    pyemu.os_utils.run("python mf6cts.py config_file.py", cwd=t_d)
    invest_run(t_d)

    pf = pyemu.utils.PstFrom(t_d, "template", spatial_reference=sim.get_model("gwt").modelgrid, remove_existing=True)
    # pf.mod_sys_cmds.append("python mf6cts.py config_file.py")
    pf.mod_py_cmds.append("run()")
    pf.add_py_function("cts_mf6_test.py", "run()", is_pre_cmd=None)
    pf.extra_py_imports.append("flopy")
    pf.extra_py_imports.append("shutil")

    dxdy = np.zeros(nrowncol) + delrdelc

    gs = pyemu.geostats.GeoStruct(variograms=pyemu.geostats.ExpVario(1.0, delrdelc * 5))
    pf.add_parameters(os.path.join(sim_ws, "gwf.npf_k.txt"), par_type="grid", par_style="d", lower_bound=0.1,
                      upper_bound=100,
                      geostruct=gs, par_name_base="hk")
    wfiles = [os.path.join(sim_ws, f) for f in os.listdir(sim_ws) if "wel_stress_period" in f and f.endswith(".txt")]
    pf.add_parameters(wfiles, par_type="constant", par_style="m", transform="none", upper_bound=2, lower_bound=0,
                      par_name_base="wmult",
                      index_cols=[0, 1, 2], use_cols=[3])

    for wfile in wfiles:
        sp = wfile.split('.')[0].split('_')[-1]
        pf.add_parameters(wfile, par_type="constant", par_style="m", transform="none", upper_bound=2, lower_bound=0,
                          pargp="wmult_sp:{0}".format(sp),
                          index_cols=[0, 1, 2], use_cols=[3])

    pf.add_observations("heads.dat", prefix="heads")
    pf.add_observations("concen.dat", prefix="concen")
    pf.add_observations("apiheads.dat", prefix="apiheads")
    pf.add_observations("apiconcen.dat", prefix="apiconcen")

    sum_files = ["gwf_cts_flow_node_summary.csv", "gwf_cts_flow_system_summary.csv",
                 "gwt_cts_node_summary.csv", "gwt_cts_system_summary.csv"]
    pf.add_observations("gwf_cts_flow_node_summary.csv", index_cols=[0, 2, 4, 7, 8], use_cols=[9, 10, 11, 12],
                        ofile_sep=",",
                        prefix="nodeflow")
    pf.add_observations("gwf_cts_flow_system_summary.csv", index_cols=[0, 2, 4], use_cols=[7, 8, 9, 10], ofile_sep=",",
                        prefix="sysflow")
    pf.add_observations("gwt_cts_node_summary.csv", index_cols=[0, 2, 4, 7], use_cols=[8, 9, 10, 11, 12],
                        ofile_sep=",",
                        prefix="nodemass")
    pf.add_observations("gwt_cts_system_summary.csv", index_cols=[0, 2, 4], use_cols=[7, 8, 9, 10, 11, 12, 13, 14, 15],
                        ofile_sep=",",
                        prefix="sysmass")

    lstcsv_files = [f for f in os.listdir("template") if "inc" in f or "cum" in f]
    for f in lstcsv_files:
        df = pd.read_csv(os.path.join("template", f))
        pf.add_observations(f, index_cols=df.columns[0], use_cols=df.columns.tolist()[1:], prefix=f.split(".")[0])

    pst = pf.build_pst()
    pe = pf.draw(num_reals=100, use_specsim=True)
    pe.enforce()
    pe.to_csv(os.path.join("template", "prior.csv"))
    pst.pestpp_options["ies_par_en"] = "prior.csv"

    pst.control_data.noptmax = 0
    pst.write(os.path.join("template", "pest.pst"))
    pyemu.os_utils.run("pestpp-ies pest.pst", cwd="template")

    pst.control_data.noptmax = -1
    pst.write(os.path.join("template", "pest.pst"))
    pyemu.os_utils.start_workers("template", "pestpp-ies", "pest.pst", num_workers=10, worker_root=".",
                                 master_dir="prior_mc_master")


def plot_prior_mc():
    """plot the prior mc results
    """
    import matplotlib.pyplot as plt
    m_d = "prior_mc_master"
    wel_data = {}

    with open(os.path.join("fivespot_api", "gwf.wel"), 'r') as f:
        while True:
            line = f.readline()
            if line == "":
                break
            if line.lower().strip().startswith("begin period"):
                wd = []
                period = int(line.strip().split()[-1])
                while True:
                    line2 = f.readline()
                    if line2 == "":
                        raise Exception()
                    if line2.lower().strip().startswith("end period"):
                        break
                    raw = line2.strip().split()
                    fname = raw[1].replace('\'',"")
                    wd = pd.read_csv(os.path.join("fivespot_api",fname),header=None,delim_whitespace=True)
                    #l, r, c = int(raw[0]) - 1, int(raw[1]) - 1, int(raw[2]) - 1
                    #flx = float(raw[3])
                    #wd.append([l, r, c, flx])
                wel_data[period] = wd.values
    wel_mask = np.zeros((nrowncol, nrowncol)) - 1
    for i,j,flx in zip(wd.iloc[:,1],wd.iloc[:,2],wd.iloc[:,3]):

        if flx > 0:
            wel_mask[i-1, j-1] = 1
        else:
            wel_mask[i-1, j-1] = 0
    wel_mask = np.ma.masked_where(wel_mask == -1, wel_mask)

    pst = pyemu.Pst(os.path.join(m_d, "pest.pst"))
    oe = pd.read_csv(os.path.join(m_d, "pest.0.obs.csv"), index_col=0)
    pe = pd.read_csv(os.path.join(m_d, "pest.0.par.csv"), index_col=0)
    par = pst.parameter_data
    hkpar = par.loc[par.parnme.str.contains("hk"), :].copy()
    hkpar.loc[:, 'i'] = hkpar.i.apply(np.int)
    hkpar.loc[:, 'j'] = hkpar.j.apply(np.int)
    hkmn = np.log10(hkpar.parlbnd.min())
    hkmx = np.log10(hkpar.parubnd.max())

    # plot a time series of mass removed
    obs = pst.observation_data
    obs.loc[:, "totim"] = obs.totim.apply(np.float)

    wobs = obs.loc[obs.obsnme.apply(lambda x: "gwtcum" in x and "wel-out" in x and "api" not in x), :].copy()
    wobs.sort_values(by="totim", inplace=True)
    awobs = obs.loc[obs.obsnme.apply(lambda x: "gwtcum" in x and "wel-out" in x and "api" in x), :].copy()
    awobs.sort_values(by="totim", inplace=True)
    tmax = wobs.totim.max()
    o1 = wobs.loc[wobs.totim == tmax, "obsnme"]
    o2 = awobs.loc[awobs.totim == tmax, "obsnme"]

    d = oe.loc[:, o1].values - oe.loc[:, o2].values
    reals = ["base", oe.index[np.argmin(d)]]

    cobs = obs.loc[obs.obsnme.apply(lambda x: "concen" in x and "arr" in x and not "api" in x), :].copy()
    cobs.loc[:, 'i'] = cobs.i.apply(np.int)
    cobs.loc[:, 'j'] = cobs.j.apply(np.int)
    hobs = obs.loc[obs.obsnme.apply(lambda x: "head" in x and "arr" in x and not "api" in x), :].copy()
    hobs.loc[:, 'i'] = hobs.i.apply(np.int)
    hobs.loc[:, 'j'] = hobs.j.apply(np.int)

    acobs = obs.loc[obs.obsnme.apply(lambda x: "concen" in x and "arr" in x and "api" in x), :].copy()
    acobs.loc[:, 'i'] = acobs.i.apply(np.int)
    acobs.loc[:, 'j'] = acobs.j.apply(np.int)
    ahobs = obs.loc[obs.obsnme.apply(lambda x: "head" in x and "arr" in x and "api" in x), :].copy()
    ahobs.loc[:, 'i'] = ahobs.i.apply(np.int)
    ahobs.loc[:, 'j'] = ahobs.j.apply(np.int)
    nrow = cobs.i.max() + 1
    ncol = cobs.j.max() + 1
    fig, axes = plt.subplots(len(reals), 3, figsize=(7, 3.7))
    # row_labs = ["{0:5.2f} mg difference".format(d.min()),"homogenous","{0:5.2f} mg difference".format(d.max())]
    row_labs = ["homogeneous", "heterogenous"]
    acount = 0

    for ireal, real_idx in enumerate(reals):
        carr = np.zeros((nrow, ncol))
        carr[cobs.i, cobs.j] = oe.loc[real_idx, cobs.obsnme]
        harr = np.zeros((nrow, ncol))
        harr[hobs.i, hobs.j] = oe.loc[real_idx, hobs.obsnme]
        acarr = np.zeros((nrow, ncol))
        acarr[acobs.i, acobs.j] = oe.loc[real_idx, acobs.obsnme]
        aharr = np.zeros((nrow, ncol))
        aharr[ahobs.i, ahobs.j] = oe.loc[real_idx, ahobs.obsnme]
        hkarr = np.zeros((nrow, ncol))
        hkarr[hkpar.i, hkpar.j] = pe.loc[real_idx, hkpar.parnme].apply(np.log10)
        c = axes[ireal, 0].imshow(hkarr, vmin=hkmn, vmax=hkmx, alpha=0.75, cmap="cividis")
        plt.colorbar(c, ax=axes[ireal, 0], label="HK ($log_{10}\\frac{m}{d}$)")
        c = axes[ireal, 1].imshow(carr, vmin=0, vmax=100, alpha=0.85, cmap="cividis")
        plt.colorbar(c, ax=axes[ireal, 1], label="concentration ($\\frac{mg}{l}$)")
        c = axes[ireal, 2].imshow(acarr, vmin=0, vmax=100, alpha=0.85, cmap="cividis")
        plt.colorbar(c, ax=axes[ireal, 2], label="concentration ($\\frac{mg}{l}$)")

        axes[ireal, 0].imshow(wel_mask, cmap="cool", vmin=0, vmax=1)
        axes[ireal, 1].imshow(wel_mask, cmap="cool", vmin=0, vmax=1)
        axes[ireal, 2].imshow(wel_mask, cmap="cool", vmin=0, vmax=1)
        cs = axes[ireal, 1].contour(harr, levels=4, colors="k", linewidths=0.25)
        plt.clabel(cs)
        cs = axes[ireal, 2].contour(aharr, levels=4, colors="k", linewidths=0.25)
        plt.clabel(cs)
        axes[ireal, 0].set_title("{0}) HK {1}".format(string.ascii_uppercase[acount], row_labs[ireal]), loc="left")
        acount += 1
        axes[ireal, 1].set_title("{0}) unbalanced flows\n{1}".format(string.ascii_uppercase[acount], row_labs[ireal]),
                                 loc="left")
        acount += 1
        axes[ireal, 2].set_title("{0}) balanced flows \n{1}".format(string.ascii_uppercase[acount], row_labs[ireal]),
                                 loc="left")
        acount += 1

    for ax in axes.flatten():
        ax.set_xlabel("column")
        ax.set_ylabel("row")
        lm = ax.get_xlim()
        ax.set_xlim(lm[0] - 2, lm[1] + 2)
        lm = ax.get_ylim()
        ax.set_ylim(lm[0] + 2, lm[1] - 2)

    plt.tight_layout()
    plt.savefig("real_summary.pdf")

    gobs = obs.loc[obs.obsnme.apply(lambda x: "gwtcum" in x and "wel-out" in x and "api" not in x), :].copy()
    gobs.sort_values(by="totim", inplace=True)
    agobs = obs.loc[obs.obsnme.apply(lambda x: "gwtcum" in x and "wel-out" in x and "api" in x), :].copy()
    agobs.sort_values(by="totim", inplace=True)

    # gobs = gobs.loc[gobs.totim==gobs.totim.max(),:]
    # agobs = agobs.loc[agobs.totim == agobs.totim.max(), :]

    goobs = obs.loc[obs.obsnme.apply(lambda x: "gwfinc" in x and "wel-out" in x and "api" not in x), :].copy()
    goobs.sort_values(by="totim", inplace=True)
    agoobs = obs.loc[obs.obsnme.apply(lambda x: "gwfinc" in x and "wel-out" in x and "api" in x), :].copy()
    agoobs.sort_values(by="totim", inplace=True)

    goobsin = obs.loc[obs.obsnme.apply(lambda x: "gwfinc" in x and "wel-in" in x and "api" not in x), :].copy()
    goobsin.sort_values(by="totim", inplace=True)
    agoobsin = obs.loc[obs.obsnme.apply(lambda x: "gwfinc" in x and "wel-in" in x and "api" in x), :].copy()
    agoobsin.sort_values(by="totim", inplace=True)

    # print(wobs)

    fig, axes = plt.subplots(2, 2, figsize=(6, 4))

    twin_ax = []
    for ireal, real in enumerate(reals):
        axes[0, ireal].plot(agoobs.totim.values, oe.loc[real, agoobs.obsnme].values, "0.5", label="extraction flux")
        axes[0, ireal].plot(agoobsin.totim.values, oe.loc[real, agoobsin.obsnme].values, "k--", label="injection flux")
        axes[1, ireal].plot(goobs.totim.values, oe.loc[real, goobs.obsnme].values, "0.5", label="extraction flux")
        axes[1, ireal].plot(goobsin.totim.values, oe.loc[real, goobsin.obsnme].values, "k--", label="injection flux")
        axt = plt.twinx(axes[0, ireal])
        axt.scatter(agobs.totim.values, oe.loc[real, agobs.obsnme].values, marker=".", c='m', label="mass extracted")
        twin_ax.append(axt)
        axt = plt.twinx(axes[1, ireal])
        axt.scatter(gobs.totim.values, oe.loc[real, gobs.obsnme].values, marker=".", c='m', label="mass extracted")
        twin_ax.append(axt)

    ymn = min([ax.get_ylim()[0] for ax in twin_ax])
    ymx = max([ax.get_ylim()[1] for ax in twin_ax])
    for ax in twin_ax:
        ax.set_ylim(ymn, ymx)
        ax.set_ylabel("cumulative mass extracted ($mg$)")
        xhand, xlab = ax.get_legend_handles_labels()

    ymn = min([ax.get_ylim()[0] for ax in axes.flatten()])
    ymx = max([ax.get_ylim()[1] for ax in axes.flatten()])
    labels = ["A) balanced flows homogenous", "B) balanced flows hetergenous",
              "C) unbalanced flows homogenous", "D) unbalanced flows heterogenous"]
    for lab, ax in zip(labels, axes.flatten()):
        ax.set_ylim(ymn, ymx)
        h, l = ax.get_legend_handles_labels()
        h.extend(xhand)
        l.extend(xlab)
        ax.legend(h, l, loc="upper left")
        ax.set_title(lab, loc="left")
        ax.set_xlabel("time (d)")
        ax.set_ylabel("WEL flux ($\\frac{m^3}{d}$)")

    plt.tight_layout()
    plt.savefig("treatment_summary.pdf")
    plt.close(fig)


def create_mt3dusgs(org_sim_ws, org_simt_ws):

    """ this never worked...

    :param org_sim_ws:
    :param org_simt_ws:
    :return:
    """
    sim = flopy.mf6.MFSimulation.load(sim_ws=org_sim_ws)
    gwf = sim.get_model("gwf")
    simt_ws = org_simt_ws + "_mt"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)

    dz = gwf.modelgrid.delz
    mt = flopy.mt3d.Mt3dms(model_ws=simt_ws, ftlfilename=None)
    btm = flopy.mt3d.Mt3dBtn(mt, nlay=gwf.dis.nlay.data, nrow=gwf.dis.nrow.data, ncol=gwf.dis.ncol.data,
                             nper=sim.tdis.nper.data,
                             delc=gwf.dis.delc.array, delr=gwf.dis.delr.array, perlen=perlen, htop=gwf.dis.top.array,
                             nstp=1, tsmult=1,
                             laycon=1, dz=dz, sconc=1)
    adv = flopy.mt3d.Mt3dAdv(mt, mixelm=0)
    gcg = flopy.mt3d.Mt3dGcg(mt)
    dsp = flopy.mt3d.Mt3dDsp(mt)

    wel_dict = gwf.wel.stress_period_data.data
    wel_kpers = list(wel_dict.keys())
    wel_kpers.sort()

    ghb_dict = gwf.ghb.stress_period_data.data
    ghb_kpers = list(ghb_dict.keys())
    ghb_kpers.sort()

    itype_dict = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_spd = {}
    for kper in range(sim.tdis.nper.data):
        ssm_data = []
        if kper in ghb_dict:
            for item in ghb_dict[kper]:
                # print(item)
                ssm_data.append([item[0][0], item[0][1], item[0][2], item[2], itype_dict["GHB"]])
        if kper in wel_dict:
            for item in wel_dict[kper]:
                # print(item)
                ssm_data.append([item[0][0], item[0][1], item[0][2], item[2], itype_dict["WEL"]])
        if len(ssm_data) > 0:
            ssm_spd[kper] = ssm_data
    ssm = flopy.mt3d.Mt3dSsm(mt, stress_period_data=ssm_spd, mxss=10000, )

    mt.write_input()

    # now add the 3 mf6 flow binary files needed for ft link
    flow_files = ["gwf.hds", "gwf.bud", "gwf.dis.grb"]
    flow_units = [20, 10, 21]
    for f in flow_files:
        shutil.copy2(os.path.join(org_sim_ws, f), os.path.join(simt_ws, f))

    nam_file = os.path.join(simt_ws, "mt3dtest.nam")
    lines = open(nam_file, 'r').readlines()
    with open(nam_file, 'w') as f:
        for line in lines:
            f.write(line)
        for u, fname in zip(flow_units, flow_files):
            f.write("FT6 {0} {1}\n".format(u, fname))
        # f.write("CTS  68  mt3dusgs.cts\n")
        # f.write("DATA  69  mt3dusgs.cto\n")

    # read the existing mf6 cts package
    cts_spd = {}
    with open(os.path.join(org_simt_ws, "model.cts"), 'r') as f:
        while True:
            line = f.readline()
            if line == "":
                break
            if "begin period" in line.lower():
                sp = int(line.split("period")[1].split()[0])
                if sp not in cts_spd:
                    cts_spd[sp] = {}
                    iline = 1
                cts = int(line.split("cts")[1].split()[0])
                eff = float(line.split("efficiency")[1].split()[0])
                inj, ext = [], []

                while True:
                    line2 = f.readline()

                    if line2 == "":
                        raise Exception("EOF for {0}".format(line))
                    if "end period" in line2.lower():
                        break
                    raw = line2.strip().split()
                    l, r, c = [int(i) for i in raw[-3:]]
                    if raw[2] == "out":
                        ext.append([l, r, c, iline])
                    else:
                        inj.append([l, r, c, iline])
                    iline += 1
                cts_spd[sp][cts] = [ext, inj, eff]
    print(cts_spd[2])

    with open(os.path.join(simt_ws, "mt3dusgs.cts"), 'w') as f:
        # MXCTS, ICTSOUT, MXEXT, MXINJ, MXWEL, IFORCE,ICTSPKG
        f.write("{0:9d} {1:9d} {2:9d} {3:9d} {4:9d} {5:9d} {6:9d}\n".format(2, 69, 4, 4, 10, 1, 1))
        for kper in range(len(perlen)):
            if kper + 1 not in cts_spd:
                # NCTS
                f.write("{0:9d}\n".format(0))
                continue
            # NCTS
            f.write("{0:9d}\n".format(len(cts_spd[kper + 1])))
            for cts, spd in cts_spd[kper + 1].items():
                # ICTS, NEXT, NINJ, ITRTINJ
                f.write("{0:9d} {1:9d} {2:9d} {3:9d}\n".format(cts, len(spd[0]), len(spd[1]), 1))
                for item in spd[0]:
                    # KEXT, IEXT, JEXT, IWEXT
                    f.write("{0:9d} {1:9d} {2:9d} {3:9d}\n".format(*item))
                # IOPTINJ CMCHINJ
                f.write("{0:9d} {1:9.3f}\n".format(1, -1. * (spd[2])))
                for item in spd[1]:
                    # KINJ, IINJ, JINJ, IWINJ
                    f.write("{0:9d} {1:9d} {2:9d} {3:9d}\n".format(*item))

    shutil.copy2(mt3d_bin, os.path.join(simt_ws, os.path.split(mt3d_bin)[-1]))
    # time.sleep(5)
    # pyemu.os_utils.run("./mt3dusgs mt3dtest.nam".format(os.path.split(mt3d_bin)[-1]),cwd=simt_ws)


def test_five_spotish_simple_api_off():
    """test all cts instances "off" for all stress periods

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.dll"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, simple_pattern=True, eff1=0.0, eff2=0.0, nlay=1)
    cts_file = os.path.join(org_sim_ws + "_t", "model.cts")

    lines = open(cts_file, 'r').readlines()
    with open(cts_file, 'w') as f:
        for line in lines:
            if line.lower().startswith("wel"):
                continue
            f.write(line)

    t = org_sim_ws + "simple"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_sim_ws, t)
    base_sim_ws = org_sim_ws
    org_sim_ws = t

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = base_sim_ws + "_t"
    t = base_sim_ws + "simple_t"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_simt_ws, t)
    org_simt_ws = t
    t = None

    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()

    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    assert np.abs(cum.loc[cum.index[-1], "wel"]) == np.abs(api_cum.loc[api_cum.index[-1], "wel"])
    abs_frac_diff = np.abs(api_cum.loc[api_cum.index[-1], "wel"] / cum.loc[cum.index[-1], "wel"])
    assert abs_frac_diff == 1.0

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    assert node_df.shape[0] == 0

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    assert sys_df.shape[0] == 0


def test_five_spotish_simple_api_off2():
    """Test cts "off" for a few stress periods

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.dll"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, simple_pattern=True, eff1=0.0, eff2=0.0, nlay=1)
    cts_file = os.path.join(org_sim_ws + "_t", "model.cts")

    lines = open(cts_file, 'r').readlines()
    period = 0
    with open(cts_file, 'w') as f:
        for line in lines:
            if line.lower().startswith("begin period"):
                period = int(line.strip().split()[2])
            if line.lower().startswith("wel") and period == 6:
                continue
            f.write(line)

    t = org_sim_ws + "simple"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_sim_ws, t)
    base_sim_ws = org_sim_ws
    org_sim_ws = t

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = base_sim_ws + "_t"
    t = base_sim_ws + "simple_t"
    if os.path.exists(t):
        shutil.rmtree(t)
    shutil.copytree(org_simt_ws, t)
    org_simt_ws = t
    t = None

    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()

    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    assert np.abs(cum.loc[cum.index[-1], "wel"]) > np.abs(api_cum.loc[api_cum.index[-1], "wel"])
    abs_frac_diff = np.abs(api_cum.loc[api_cum.index[-1], "wel"] / cum.loc[cum.index[-1], "wel"])

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    assert 6 not in node_df.stress_period.astype(int).values

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    assert 6 not in sys_df.stress_period.astype(int).values


def test_five_spotish_api_shared_inj():
    """a basic test of the cts

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)

    setup_five_spotish(plot=False, sim_ws=org_sim_ws, eff1=0.7,eff2=0.7, simple_eff=True, nlay=1,
                       ghb_source=100,
                       simple_pattern=True, simple_hk=True,shared_inj=True)

    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = org_sim_ws + "_t"
    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    flow_budget_file = "gwf.bud"
    # shutil.copy2(os.path.join(sim_ws,flow_budget_file),os.path.join(simt_ws,flow_budget_file))
    try:
        mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                     is_structured=True)
    except Exception:
        return
    else:
        raise Exception("should have failed")


    mf6.solve_gwf()
    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0, axis=1),
                   :].mass.sum()
    out_node_mass = node_df.loc[node_df.apply(lambda x: x.cum_vol < 0, axis=1),
                    :].mass.sum()
    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    if in_node_mass != 0.0:
        assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    eff = sys_df.mass_treated / sys_df.mass_removed
    d = (eff - sys_df.requested_efficiency).apply(np.abs)
    print(d)
    assert d.max() < 1.0e-10

    d = (sys_df.mass_removed - (sys_df.mass_treated + sys_df.mass_injected)).apply(np.abs)
    print(d)

    wel_df = pd.read_csv(os.path.join(sim_ws,"gwf.wel_stress_period_data_2.txt"),header=None,delim_whitespace=True,names=["l","r","c","flux","concen"])
    print(wel_df)

    node_df = pd.read_csv(os.path.join(sim_ws, "gwf_cts_flow_node_summary.csv"))
    cts_vals = node_df.loc[node_df.stress_period==2,"requested_rate"]
    cts_vals.values.sort()
    wel_vals = wel_df.flux.values
    wel_vals.sort()
    d = np.abs(wel_vals - cts_vals).sum()
    print(d)
    assert d < 1.0e-10


def test_five_spotish_api_maw_configfile_staggered():
    """Test MAW with the config file usage by solving flow then transport in steps

    """

    from mf6cts import Mf6Cts

    # the mf6 library

    # lib_name = "libmf6.so"
    # lib_path = os.path.join(".", lib_name)
    # the model files directory
    org_sim_ws = "fivespot"
    np.random.seed(111)
    setup_five_spotish(plot=False, sim_ws=org_sim_ws, nlay=1)
    convert_5spot_to_maw("fivespot", nlay=1)
    org_sim_ws = org_sim_ws + "_maw"
    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    # this dir should have been created with the call to setup_five_spotish()
    org_simt_ws = org_sim_ws + "_t"
    assert os.path.exists(org_simt_ws), org_simt_ws

    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    # first copy the mf6cts.py file to this dir
    src_file = os.path.join("..", "mf6cts", "mf6cts.py")
    dest_file = "mf6cts.py"
    if os.path.exists(dest_file):
        os.remove(dest_file)
    shutil.copy2(src_file, dest_file)

    # now write a config file
    config_file = "config_file.py"
    if os.path.exists(config_file):
        os.remove(config_file)
    with open(config_file, 'w') as f:
        f.write("cts_filename='{0}'\n".format("model.cts"))
        f.write("lib_name='{0}'\n".format(os.path.split(lib_name)[-1]))
        f.write("transport_dir='{0}'\n".format(simt_ws))
        f.write("flow_dir='{0}'\n".format(sim_ws))
        f.write("is_structured=True\n")
        f.write("flow_output_files=['gwf.hds','gwf.bud','gwf.maw.bud']\n")
        f.write("solve_gwf=True\n")
        f.write("transfer_flow_output_files=True\n")
        f.write("solve_gwt=False\n")

    for fname in ['gwf.hds', 'gwf.bud', 'gwf.maw.bud']:
        if os.path.exists(os.path.join(simt_ws, fname)):
            os.remove(os.path.join(simt_ws, fname))

    os.system("python mf6cts.py config_file.py")
    #os.remove(dest_file)

    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    # assert np.abs(cum.loc[cum.index[-1], "maw"]) > np.abs(api_cum.loc[api_cum.index[-1], "maw"])
    abs_fr_diff = np.abs(api_cum.loc[api_cum.index[-1], "maw"]) / np.abs(cum.loc[cum.index[-1], "maw"])
    assert abs_fr_diff < 0.025, abs_fr_diff

    node_df = pd.read_csv(os.path.join(sim_ws, "gwf_cts_flow_node_summary.csv"))
    node_df = node_df.loc[node_df.requested_rate > 0, :]
    d = (node_df.requested_rate - node_df.actual_rate).apply(np.abs).sum()
    print(d)
    assert d > 1000



    # make sure the gwt stuff is empty
    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    assert node_df.shape[0] == 0

    # now write a config file for only solving gwt
    config_file = "config_file.py"
    if os.path.exists(config_file):
        os.remove(config_file)
    with open(config_file, 'w') as f:
        f.write("cts_filename='{0}'\n".format("model.cts"))
        f.write("lib_name='{0}'\n".format(os.path.split(lib_name)[-1]))
        f.write("transport_dir='{0}'\n".format(simt_ws))
        f.write("flow_dir='{0}'\n".format(sim_ws))
        f.write("is_structured=True\n")
        f.write("flow_output_files=['gwf.hds','gwf.bud','gwf.maw.bud']\n")
        f.write("solve_gwf=False\n")
        f.write("transfer_flow_output_files=False\n")
        f.write("solve_gwt=True\n")

    os.system("python mf6cts.py config_file.py")
    os.remove(dest_file)

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01, abs_frac_diff

    abs_frac_diff = np.abs((in_node_mass - api_cum.MWT_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.MWT_IN.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

def fr1_test():
    from mf6cts import Mf6Cts
    # instaniate...
    sim_ws = "fr1_test"
    sim = flopy.mf6.MFSimulation(sim_name="mfsim", sim_ws=sim_ws, continue_=True, memory_print_option="all")
    perlen = [10000.0 for _ in range(3)]

    nlay,nrow,ncol = 1,1,5
    top = 1
    botm = 0
    # build up the time stepping container for the tdis package
    tdis_data = [(p, 1, 1.0) for p in perlen]

    tdis = flopy.mf6.ModflowTdis(simulation=sim, nper=len(tdis_data), perioddata=tdis_data)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions="newton")

    # instantiate discretization package
    id = np.ones((nrow,ncol),dtype=int)
    id[0,2:-2] = 0


    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=1, delc=1,
                                  top=top,
                                  botm=botm,idomain=id)

    # instantiate node property flow package
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=100, icelltype=1,
                                  save_specific_discharge=True,
                                  save_flows=True, save_saturation=True)

    # instantiate initial conditions for the flow solution - set starting heads at midpoint of layer 1
    ic = flopy.mf6.ModflowGwfic(gwf, strt=top)

    # instantiate the storage package - stress period 1 is steady state, transient after that...
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1, steady_state={0: True}, transient={1: True}, ss=0.00001, sy=0.01)

    # output control - headsave and budget file names
    oc = flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=bud_file, head_filerecord=hds_file,
                                headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
                                saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
                                printrecord=[("BUDGET", "LAST"),("HEAD","LAST")], )


    wel_data = [[(0, 0, 1), -10.0, 0.0]]
    wel_data.append([(0,0,ncol-2),10.0,0.0])
    # instantiate the wel package
    wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0:wel_data}, save_flows=True,
                                  auxiliary="concentration", auto_flow_reduce=1.0)


    ghb_data = []
    ghb_data.extend([[(0, 0, 0), top, 1000.0, 1.0]])
    ghb_data.extend([[(0, 0,ncol-1), (top + botm) / 2., 1000., 0.0]])
    #ghb_data.extend([[(0, 0, int(ncol/2)), top, 100., 0.0]])

    ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghb_data, auxiliary="concentration",
                                  save_flows=True)

    # just a generic solver will do
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab", outer_dvclose=0.0001, inner_dvclose=0.0001,
                               outer_maximum=100,
                               inner_maximum=250)

    # write the input file
    sim.simulation_data.max_columns_of_data = sim.get_model("gwf").dis.nrow.data
    sim.set_all_data_external(check_data=False)
    sim.write_simulation()

    # run the flow model
    # os.chdir(sim_ws)
    shutil.copy2(mf6_bin, os.path.join(sim_ws, os.path.split(mf6_bin)[1]))
    pyemu.os_utils.run("mf6", cwd=sim_ws)
    # os.chdir("..")

    # now make a separate transport simulation and model
    simt_ws = sim_ws + "_t"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    os.makedirs(simt_ws)
    # transport sim
    simt = flopy.mf6.MFSimulation(sim_ws=simt_ws, memory_print_option="all", continue_=True)
    # transport sim temporal discet matches flow solution discret
    tdist = flopy.mf6.ModflowTdis(simulation=simt, nper=len(tdis_data), perioddata=tdis_data)

    # transport model instance
    gwtname = "gwt"
    gwt = flopy.mf6.ModflowGwt(simt, modelname=gwtname, save_flows=True)

    # transport model discret package
    dist = flopy.mf6.ModflowGwtdis(gwt, nlay=nlay, nrow=nrow, ncol=ncol, delr=1, delc=1,
                                   top=top,
                                   botm=botm,idomain=id)

    # copy in the flow solution output files for use in the transport model
    shutil.copy2(os.path.join(sim_ws, hds_file), os.path.join(simt_ws, hds_file))
    shutil.copy2(os.path.join(sim_ws, bud_file), os.path.join(simt_ws, bud_file))
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=[["gwfhead", hds_file], ["gwfbudget", bud_file]],
                                  flow_imbalance_correction=True)

    # initial concen
    strt = 0.0
    ict = flopy.mf6.ModflowGwtic(gwt, strt=strt)

    # remaining transport packages
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.0001)
    dsp = flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, alh=1.0, ath1=0.1, ath2=0.1)
    ssm = flopy.mf6.ModflowGwtssm(gwt,
                                  sources=[["WEL_0", "AUX", "CONCENTRATION"], ["GHB_0", "AUX", "CONCENTRATION"]])

    # transport sim output
    ucn_file = "{}.ucn".format(gwt.name)
    oct = flopy.mf6.ModflowGwtoc(gwt, budget_filerecord="{}.cbc".format(gwt.name),
                                 concentration_filerecord=ucn_file,
                                 concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
                                 saverecord=[("CONCENTRATION", "ALL")],
                                 printrecord=[("BUDGET", "LAST"),("CONCENTRATION","LAST")])

    imst = flopy.mf6.ModflowIms(simt, filename="{}.ims".format(gwt.name), linear_acceleration="bicgstab",
                                inner_dvclose=0.0001, outer_dvclose=0.0001, outer_maximum=100, inner_maximum=100)

    simt.register_ims_package(imst, [gwt.name])

    # write a cts file
    effs = [0.0,0.5,1.0]
    with open(os.path.join(simt_ws, "model.cts"), 'w') as f:
        f.write("begin options\n\nend options\n\n")
        for kper,eff in enumerate(effs):
            f.write("begin period {0} cts 1 efficiency {1:4.3f}\n".format(kper+1, eff))
            for wd in wel_data:
                f.write(
                    "wel wel_0 {0} {1} {2} {3}\n".format("out" if wd[1] < 0 else "in", wd[0][0] + 1, wd[0][1] + 1,
                                                         wd[0][2] + 1))
            f.write("end period {0} cts 1\n\n".format(kper + 2))

    # write the transport inputs and run
    simt.write_simulation()
    shutil.copy2(mf6_bin, os.path.join(simt_ws, os.path.split(mf6_bin)[1]))
    pyemu.os_utils.run("mf6", cwd=simt_ws)

    # the api part

    org_sim_ws = sim_ws
    sim_ws = org_sim_ws + "_api"
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(org_sim_ws, sim_ws)
    shutil.copy2(lib_name, os.path.join(sim_ws, os.path.split(lib_name)[-1]))

    org_simt_ws = simt_ws
    simt_ws = org_simt_ws + "_api"
    if os.path.exists(simt_ws):
        shutil.rmtree(simt_ws)
    shutil.copytree(org_simt_ws, simt_ws)
    shutil.copy2(lib_name, os.path.join(simt_ws, os.path.split(lib_name)[-1]))

    mf6 = Mf6Cts("model.cts", os.path.split(lib_name)[-1], transport_dir=simt_ws, flow_dir=sim_ws,
                 is_structured=True)

    mf6.solve_gwf()

    shutil.copy2(os.path.join(sim_ws, "gwf.hds"), os.path.join(simt_ws, "gwf.hds"))
    shutil.copy2(os.path.join(sim_ws, "gwf.bud"), os.path.join(simt_ws, "gwf.bud"))
    mf6.solve_gwt()
    mf6.finalize()

    # first copy the mf6cts.py file to this dir
    src_file = os.path.join("..", "mf6cts", "mf6cts.py")
    dest_file = "mf6cts.py"
    if os.path.exists(dest_file):
        os.remove(dest_file)
    shutil.copy2(src_file, dest_file)

    # now write a config file
    config_file = "config_file.py"
    if os.path.exists(config_file):
        os.remove(config_file)
    with open(config_file, 'w') as f:
        f.write("cts_filename='{0}'\n".format("model.cts"))
        f.write("lib_name='{0}'\n".format(os.path.split(lib_name)[-1]))
        f.write("transport_dir='{0}'\n".format(simt_ws))
        f.write("flow_dir='{0}'\n".format(sim_ws))
        f.write("is_structured=True\n")
        f.write("flow_output_files=['gwf.hds','gwf.bud']\n")

    for fname in ['gwf.hds', 'gwf.bud', 'gwf.maw.bud']:
        if os.path.exists(os.path.join(simt_ws, fname)):
            os.remove(os.path.join(simt_ws, fname))

    os.system("python mf6cts.py config_file.py")
    os.remove(dest_file)


    mf6 = None
    api_lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, gwfname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=True)
    lst = flopy.utils.Mf6ListBudget(os.path.join(org_sim_ws, gwfname + ".lst"))
    inc, cum = lst.get_dataframes(diff=True)

    print(api_cum.iloc[-1, :])
    print(cum.iloc[-1, :])
    assert np.abs(cum.loc[cum.index[-1], "wel"]) > np.abs(api_cum.loc[api_cum.index[-1], "wel"])
    abs_frac_diff = np.abs(api_cum.loc[api_cum.index[-1], "wel"] / cum.loc[cum.index[-1], "wel"])
    assert abs_frac_diff < 0.01

    api_lst = Mf6TListBudget(os.path.join(simt_ws, gwtname + ".lst"))
    api_inc, api_cum = api_lst.get_dataframes(diff=False)
    print(api_cum.iloc[-1, :])

    node_df = pd.read_csv(os.path.join(sim_ws, "gwf_cts_flow_node_summary.csv"))
    assert node_df.shape[0] == 6

    node_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_node_summary.csv"))
    assert node_df.shape[0] == 6
    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == 1, axis=1),
                   :].cum_mass.sum()

    out_node_mass = node_df.loc[
                    node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == 1, axis=1),
                    :].cum_mass.sum()
    print(in_node_mass, out_node_mass)
    abs_frac_diff = np.abs(in_node_mass - out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    in_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()
    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    out_node_mass = node_df.loc[
                   node_df.apply(lambda x: x.cum_vol < 0 and x.stress_period == node_df.stress_period.max(), axis=1),
                   :].cum_mass.sum()
    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    sys_df = pd.read_csv(os.path.join(simt_ws, "gwt_cts_system_summary.csv"))
    in_node_mass = sys_df.loc[
                   sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                   :].cum_mass_injected.sum()

    out_node_mass = sys_df.loc[
                    sys_df.apply(lambda x: x.cum_vol > 0 and x.stress_period == sys_df.stress_period.max(), axis=1),
                    :].cum_mass_removed.sum()

    abs_frac_diff = np.abs((in_node_mass - out_node_mass) / out_node_mass)
    print(abs_frac_diff)
    
    abs_frac_diff = np.abs((in_node_mass - api_cum.WEL_IN.max()) / in_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    abs_frac_diff = np.abs((out_node_mass - api_cum.WEL_OUT.max()) / out_node_mass)
    print(abs_frac_diff)
    assert abs_frac_diff < 0.01

    ucn = flopy.utils.HeadFile(os.path.join(simt_ws,"gwt.ucn"),text="concentration")
    for ikper,eff in enumerate(effs):
        arr = ucn.get_data(kstpkper=(0,ikper))
        sim_eff = arr[0,0,1] - (arr[0,0,-2]/arr[0,0,1])
        d = np.abs(sim_eff - eff)
        #print(arr[0,0,1],arr[0,0,-2],eff,sim_eff)
        print(d)
        assert d < 1.0e-6



if __name__ == "__main__":
    fr1_test()
    #test_five_spotish_simple_api1()
    # test_five_spotish_simple_api2()
    # plot("fivespotsimple_t", "fivespotsimple")
    # plot("fivespotsimple_t_api", "fivespotsimple_api")

    #test_five_spotish_api()

    # plot("fivespot_t", "fivespot")
    # plot("fivespot_t_api", "fivespot_api")
    # test_five_spotish_simple_api_mk2k_compare()
    #test_five_spotish_complex_api_mk2k_compare()
    # plot("fivespotsimple_t", "fivespotsimple")
    # plot("fivespotsimple_t_api", "fivespotsimple_api")
    # test_five_spotish_api_maw()
    #test_five_spotish_api_maw_complex()
    # plot("fivespot_maw_t", "fivespot_maw")
    # plot("fivespot_maw_t_api", "fivespot_maw_api")
    #test_five_spotish_api_maw_configfile_staggered()
    #plot_domain("fivespot_api")
    # plot_results_pub("fivespot_t_api", "fivespot_api")
    #setup_pst()
    #plot_prior_mc()

    # create_mt3dusgs("fivespot_api","fivespot_t_api")
    # pyemu.os_utils.run("./mt3dusgs mt3dtest.nam".format(os.path.split(mt3d_bin)[-1]), cwd="fivespot_t_api_mt")
    # test_five_spotish_simple_api_off()
    # test_five_spotish_simple_api_off2()
    #test_five_spotish_api()
    #test_five_spotish_api_shared_inj()

