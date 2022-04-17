import os
import sys
import shutil
import warnings
from datetime import datetime
import numpy as np
import modflowapi
import flopy

DT_FMT = "%Y-%m-%d %H:%M:%S"

"""
============
example CTS input file using wel:
============

begin options

end options

begin period 2 cts 1 efficiency 0.75
wel wel_0 out 3 12 12
wel wel_0 in 1 1 1
wel wel_0 in 1 33 1
end period 2 cts 1

begin period 2 cts 2 efficiency 0.5
wel wel_0 out 3 17 17
wel wel_0 in 1 33 33
wel wel_0 in 1 1 33
end period 2 cts 2


============
example CTS input file using maw:
============

begin options

end options

begin period 2 cts 1 efficiency 0.5
  maw MAW_0 in 1
  maw MAW_0 out 5
end period 2

begin period 2 cts 2 efficiency 0.1
  maw MAW_0 in 2
  maw MAW_0 in 3
  maw MAW_0 in 4
  maw MAW_0 out 6
end period 2

============
example config file for command line usage (a basic python source file):
============

cts_filename='model.cts'
lib_name='libmf6.so'
transport_dir='fivespot_maw_t_api'
flow_dir='fivespot_maw_api'
is_structured=True
flow_output_files=['gwf.hds','gwf.bud','gwf.maw.bud']



"""


class CtsRecord(object):
    """simple container class to hold a single CTS record

    Parameters:
        index (int): an index of this records location in a MODFLOW-6 vector.
            For non-maw entries, this is zero-based node number.  For maw, this is
            one-based "wellno"
        inout (str): string to identify if this record is for an injector ("in") or
            extractor ("out")
        pakname (str): the standard MODFLOW-6 package name (e.g. "WEL","MAW","DRN", etc)
        instance (str): the package instance name (e.g. "WEL_0", "MAW-1", etc)


    """

    def __init__(self, index, inout, pakname, instance, extra=""):
        self.index = int(index)
        self.inout = inout.lower()
        self.pakname = pakname.lower()
        self.instance = instance.lower()
        self.extra = extra
        self.occurence = 0

    def get_tuple(self):
        return tuple([self.index,self.inout,self.pakname,self.instance])

    def get_tuple_wo_instance(self):
        return tuple([self.index,self.inout,self.pakname])

class CtsSystem(object):
    """class representing stress-period specific CTS system blocks

    Parameters:
        cts_system_number (int): the cts system number read from the CTS file
            period block
        period (int): one-based stress period number
        efficiency (float): the system efficiency
        entries (list): a list of `CtsRecord` instances


    """

    def __init__(self, cts_system_number, period, efficiency, concentration, entries):
        self.cts_system_number = cts_system_number
        self.period = period
        self._efficiencies = {}
        self._concentrations = {}
        self._entries = {}
        self.add_period_entries(period, efficiency, concentration, entries)
        self.cum_mass_ext = 0.
        self.cum_mass_treat = 0.
        self.cum_mass_inj = 0.
        self.cum_vol = 0.
        self.num_inj = 0
        self.num_ext = 0
        self.inj_rate = 0
        self.mass_removed = 0
        self.mass_treated = 0
        self.mass_inj = 0
        self.inj_conc = 0

        self.cum_req_vol = 0
        self.cum_act_vol = 0
        self.req_rate = 0
        self.act_rate = 0
        self._balance_flows = True


    def clear_metrics(self):
        """clear the incremental metrics

        """
        self.num_inj = 0
        self.num_ext = 0
        self.inj_rate = 0
        self.mass_removed = 0
        self.mass_treated = 0
        self.mass_inj = 0
        self.inj_conc = 0

        self.req_rate = 0
        self.act_rate = 0


    def add_period_entries(self, period, efficiency, concentration, entries):
        """add a new period block for this CTS system

        Parameters:
            period (int): one-based stress period number
            efficiency (float): system efficency for this period
            concentrations (float): effluent concentration for this period
            entries (list): list of `CtsRecord` instances for this period

        """
        if period in self._efficiencies:
            raise Exception("CtsSystem exception: period {0} already in eff ")
        if period in self._entries:
            raise Exception("CtsSystem exception: period {0} already in cts entries ")

        self._entries[period] = [CtsRecord(*e) for e in entries]
        self._efficiencies[period] = efficiency
        self._concentrations[period] = concentration
        self._check_for_duplicates(period)

    def _check_for_duplicates(self,period):
        cts_records = self._entries[period]
        tups = []
        for cts_rec in cts_records:
            tups.append(cts_rec.get_tuple())
        stups = set(tups)
        if len(stups) != len(tups):
            raise Exception("CtsSystem error: duplicate entries for period {0}".format(period))


    def get_eff_concen_and_entries(self, period):
        """get the efficiency and `CtsRecord` instances for a specified period
        using the MODFLOW re-use last entry convention

        Parameters:
            period (int): the one-based stress period to get efficiency and entries for.
                If `period` is less than the first period of data, then `0,[]` is returned.
                Otherwise the first period of data greater than or equal to `period` is returned

        Returns:
             int,list: the efficiency and `CtsRecord` instances corresponding to `period`

        """
        pperiod = self._find_nearest_period(period)
        if pperiod == -1:
            return 0, 1.0e+30, []
        else:
            return self._efficiencies[pperiod], self._concentrations[pperiod], self._entries[pperiod].copy()

    def _find_nearest_period(self, period):
        keys = list(self._entries.keys())
        keys.sort()
        last_key = -1
        for key in keys:
            if key > period:
                break
            last_key = key
        return last_key

    def contains_maw(self):
        for p, items in self._entries.items():
            for item in items:
                if item.pakname == "maw":
                    return True
        return False


class Mf6Cts(object):
    """class to drive both the flow and transport MODFLOW-6 models to
    implement a flow-balanced treatment system with time-varying injectors,
    extractors, and system efficiency

    Parameters:
        cts_filename (str): the CTS system filename.  This file is expected to be found in
            the `transport_dir` directory
        lib_name (str): the MODFLOW-6 shared/DLL library filename
        transport_dir (str): the transport model directory
        flow_dir (str): the flow model directory
        is_structured (bool): flag to indicate a structured grid


    """

    def __init__(self, cts_filename, lib_name, transport_dir, flow_dir, is_structured=True):
        """
        todo: setup a log file to record what is happening


        """

        # make sure the cts input file exists
        if not os.path.exists(os.path.join(transport_dir, cts_filename)):
            raise Exception("cts_filename '{0}' not found".format(os.path.join(transport_dir, cts_filename)))
        self.cts_filename = os.path.join(transport_dir, cts_filename)

        # process the flow model
        # make sure the lib exists
        if not os.path.exists(os.path.join(flow_dir, lib_name)):
            raise Exception("MODFLOW-6 shared library  '{0}' not found in flow_dir {1}".
                            format(lib_name, flow_dir))
        # find the model name
        self._gwf_model_dict, namfile_dict = Mf6Cts.get_model_names_from_mfsim(flow_dir)
        if len(self._gwf_model_dict) != 1:
            raise Exception("only one gwf model is current supported")
        self._gwf_name = list(self._gwf_model_dict.keys())[0]
        if self._gwf_model_dict[self._gwf_name] != "gwf6":
            raise Exception("model in flow_dir is not a gwf6 type: {0}". \
                            format(self._gwf_model_dict[self._gwf_name]))



        # get a list of flow term package from the flow model nam file - needed to work out the FMI
        # flow term ordering (until strings are supported in the API)
        self._gwf_ft_pak_dict = Mf6Cts.get_gwf_namfile_entry_dict(os.path.join(flow_dir, namfile_dict[self._gwf_name]))
        self._gwf = None
        self._lib_name = lib_name
        self._flow_dir = flow_dir
        self._gwf = self._initialize_gwf(lib_name,flow_dir)


        # process the transport model
        # make sure the lib exists
        if not os.path.exists(os.path.join(transport_dir, lib_name)):
            raise Exception("MODFLOW-6 shared library  '{0}' not found in transport_dir {1}".
                            format(lib_name, transport_dir))
        # find the model name
        self._gwt_model_dict, namfile_dict = Mf6Cts.get_model_names_from_mfsim(transport_dir)
        if len(self._gwt_model_dict) != 1:
            raise Exception("only one gwt model is current supported")
        self._gwt_name = list(self._gwt_model_dict.keys())[0]
        if self._gwt_model_dict[self._gwt_name] != "gwt6":
            raise Exception("model in transport_dir is not a gwt6 type: {0}". \
                            format(self._gwt_model_dict[self._gwt_name]))
        self._gwt = None
        self._transport_dir = transport_dir
        #self._gwt = self._initialize_gwt(lib_name,transport_dir)

        # containers to hold the cts info
        self._cts_instances = {}
        self._structured_mg = None
        self.is_structured = is_structured
        if self.is_structured:
            nlay = self._gwf.get_value(self._gwf.get_var_address("NLAY", self._gwf_name.upper(), "DIS"))[0]
            nrow = self._gwf.get_value(self._gwf.get_var_address("NROW", self._gwf_name.upper(), "DIS"))[0]
            ncol = self._gwf.get_value(self._gwf.get_var_address("NCOL", self._gwf_name.upper(), "DIS"))[0]
            self._structured_mg = flopy.discretization.StructuredGrid(nrow=nrow,
                                                                      ncol=ncol,
                                                                      nlay=nlay)

        else:
            raise NotImplementedError("only structured is currently supported")



        self._flow_node_summary = open(os.path.join(flow_dir, self._gwf_name + "_cts_flow_node_summary.csv"), 'w')
        self._flow_node_summary.write("stress_period,time_step,ctime,dt,cts_system,package,instance,inout,index,")
        self._flow_node_summary.write("requested_rate,actual_rate,requested_cum_vol,actual_cum_vol\n")

        self._flow_system_summary = open(os.path.join(flow_dir, self._gwf_name + "_cts_flow_system_summary.csv"), 'w')
        self._flow_system_summary.write("stress_period,time_step,ctime,dt,cts_system,num_injectors,num_extractors,")
        self._flow_system_summary.write("requested_rate,actual_rate,requested_cum_vol,actual_cum_vol\n")

        # prepare the flow summary file

        self._cts_system_summary = open(os.path.join(transport_dir, self._gwt_name + "_cts_system_summary.csv"), 'w')
        self._cts_system_summary.write("stress_period,time_step,ctime,dt,cts_system,num_injectors,num_extractors,flow_rate,cum_vol,mass_removed,")
        self._cts_system_summary.write("cum_mass_removed,mass_treated,cum_mass_treated,mass_injected,cum_mass_injected,concen_injected,requested_efficiency\n")

        self._cts_node_summary = open(os.path.join(transport_dir, self._gwt_name + "_cts_node_summary.csv"), 'w')
        self._cts_node_summary.write("stress_period,time_step,ctime,dt,cts_system,package,instance")
        self._cts_node_summary.write(",index,flow_rate,cum_vol,concen,mass_rate,mass,cum_mass\n")
        self._cts_node_record_dict = {}
        self._maw_node_record_dict = {}
        self._cts_current_nodes = []
        self._maw_current_nodes = []

        # read the cts file
        self._read_cts_input_file()

        # check if maw in a cts system
        contains_maw = False
        for cts_num, cts_system in self._cts_instances.items():
            if cts_system.contains_maw():
                contains_maw = True
                break
        self._mwt_name = None
        self._maw_node_summary = None
        # if maw is present make sure mwt is present
        if contains_maw:
            self._mwt_name = self.get_gwt_mwt_instance_name(os.path.join(transport_dir, namfile_dict[self._gwt_name]))
            if self._mwt_name is None:
                raise Exception("MAW found in cts system, but MWT package not found in transport model")
            #prep the maw-specific node csv reporting
            self._maw_node_summary = open(os.path.join(transport_dir,"gwt_maw_node_summary.csv"),'w')
            self._maw_node_summary.write("stress_period,time_step,ctime,dt,cts_system,package,instance")
            self._maw_node_summary.write(",maw_wellno,index,flow_rate,cum_vol,concen,mass_rate,mass,cum_mass\n")

        #now turn off gwt so that it releases the file handles for the
        # flow model output files - these need to be copied in to the transport dir
        # after the flow model runs
        # self._gwt.finalize()
        # self._gwt = None
        self._initial_req_rates = {}

    def _initialize_gwf(self,lib_name,flow_dir):
        # instantiate the flow model api
        gwf = modflowapi.ModflowApi(os.path.join(flow_dir, lib_name), working_directory=flow_dir)
        gwf.initialize()
        return gwf

    def _initialize_gwt(self,lib_name,transport_dir):
        # instantiate the transport model api
        gwt = modflowapi.ModflowApi(os.path.join(transport_dir, lib_name), working_directory=transport_dir)
        gwt.initialize()
        return gwt

    def finalize(self):
        """close the api and file handles

        """
        try:
            self._gwf.finalize()
        except:
            pass
        try:
            self._gwt.finalize()
        except:
            pass

        try:
            self._cts_system_summary.close()
        except:
            pass
        try:
            self._cts_node_summary.close()
        except:
            pass

        try:
            self._cts_flow_node_summary.close()
        except:
            pass

        try:
            self._cts_flow_system_summary.close()
        except:
            pass

        try:
            self._maw_node_summary.close()
        except:
            pass

    @staticmethod
    def get_gwf_namfile_entry_dict(gwf_namefile):
        """get the flow-term packages from the flow model nam file

        Parameters:
            gwf_namefile (str): the fully pathed flow model name file

        Returns:
            dict: package-order pairs (e.g. {"wel":1,"maw":2})

        """
        pak_count = 1
        gwf_ft_pak_dict = {}
        with open(gwf_namefile, 'r') as f:
            for line in f:
                if line.strip() == "":
                    continue
                if line.strip()[0] == "#":
                    continue
                if "mvr" in line.lower():
                    raise Exception("mover not supported")
                raw = line.strip().lower().split()
                if raw[0] in ["wel6", "ghb6", "maw6", "chd6", "evt6",
                              "evta6", "lak6", "rch6", "rcha6", "riv6",
                              "sfr6", "uzf6"]:
                    bud_name = raw[0].replace("6", "")
                    if bud_name in gwf_ft_pak_dict:
                        raise Exception("duplicate package name entries for {0}".format(bud_name))
                    gwf_ft_pak_dict[bud_name] = pak_count
                    pak_count += 1
        return gwf_ft_pak_dict

    @staticmethod
    def get_gwt_mwt_instance_name(gwt_namefile):
        """get the MWT package name from the transport model nam file

        Parameters:
            gwt_namefile (str): the fully pathed transport model name file

        Returns:
            str: MWT package tag

        """
        mwt_name = None
        with open(gwt_namefile, 'r') as f:
            for line in f:
                if line.strip() == "":
                    continue
                if line.strip()[0] == "#":
                    continue
                if "mwt" in line.lower():
                    raw = line.strip().lower().split()
                    if len(raw) == 2:
                        mwt_name = raw[0]
                    elif len(raw) > 2:
                        mwt_name = raw[2]
                    else:
                        raise Exception("something is wrong")

        # if mwt_name is None:
        #    raise Exception("MWT package not found in gwt_namefile '{0}'".format(gwt_namefile))
        return mwt_name

    @staticmethod
    def get_model_names_from_mfsim(sim_ws):
        """return the model names from an mfsim.nam file

        Parameters:
            sim_ws (str): the simulation path

        Returns:
            dict,dict: a pair of dicts, first is model-name:model-type (e.g. {"gwf-1":"gwf"},
                the second is model namfile: model-type (e.g. {"gwf-1":"gwf_1.nam"})

        """
        sim_nam = os.path.join(sim_ws, "mfsim.nam")
        if not os.path.exists(sim_nam):
            raise Exception("simulation nam file '{0}' not found".format(sim_nam))
        model_dict = {}
        namfile_dict = {}
        with open(sim_nam, 'r') as f:
            while True:
                line = f.readline()
                if line == "":
                    raise EOFError("EOF when looking for 'models' block")
                if line.strip().lower().startswith("begin") and "models" in line.lower():
                    while True:
                        line2 = f.readline()
                        if line2 == "":
                            raise EOFError("EOF when reading 'models' block")
                        elif line2.strip().lower().startswith("end") and "models" in line2.lower():
                            break
                        raw = line2.strip().lower().split()
                        if raw[-1] in model_dict:
                            raise Exception("duplicate model name found: '{0}'".format(raw[-1]))
                        model_dict[raw[2]] = raw[0]
                        namfile_dict[raw[2]] = raw[1]

                    break
        return model_dict, namfile_dict

    def _balance_cts_flows(self, stress_period, dt, kiter, record=False):
        """use the current simulated extraction rates for each cts system to
        balance the injection rates

        Args:
            stress_period (int): the one-based stress period number
            dt (float): time length
            kiter (int): zero-based iteration number
            record (bool): flag to store results from the is balancing

        """

        # for each cts system
        for cts_num, cts_instance in self._cts_instances.items():

            # work out what flow rates we are actually getting from the  model
            out_tot = 0.0
            req_tot = 0.0
            num_inj = 0
            num_ext = 0
            inj_nodes = []
            inj_req = 0.0
            eff, concen, cts_data = cts_instance.get_eff_concen_and_entries(stress_period)
            inj_idx = []
            cts_instance.clear_metrics()
            for cdata in cts_data:
                req_rate = None
                act_rate = None
                if cdata.pakname.lower() == "maw":

                    # requested flow per maw
                    addr = ["RATE", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    req_rates = self._gwf.get_value(wbaddr)

                    # if no flux, continue
                    if np.sum(np.abs(req_rates)) == 0:
                        continue

                    # the simulated rates from the last solve
                    addr = ["RATESIM", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    sim_rates = self._gwf.get_value(wbaddr)

                    # handle inj vs ext wells
                    req_rate = req_rates[cdata.index - 1]
                    if kiter == 1:
                        self._initial_req_rates[tuple([cts_num,cdata.index])] = req_rate
                    act_rate = sim_rates[cdata.index - 1]
                    if cdata.inout == "out":
                        req_tot += req_rate
                        out_tot += act_rate
                        num_ext += 1
                    elif cdata.inout == "in":
                        num_inj += 1
                        inj_nodes.append(cdata.index)
                        inj_req += req_rate

                else:

                    # get the requested rates - from the input file
                    addr = ["BOUND", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    req_rates = self._gwf.get_value(wbaddr)

                    # if no flux, continue
                    if np.sum(np.abs(req_rates)) == 0:
                        continue

                    # get the bc package rhs info - this has what was actually taken from the model
                    addr = ["RHS", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    sim_rates = self._gwf.get_value(wbaddr)

                    # get the node list in case things are not in order?
                    addr = ["NODELIST", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    nodes = self._gwf.get_value(wbaddr)

                    # if this is an extraction node, accum req and actual rates
                    req_rate = req_rates[np.where(nodes == cdata.index + 1)[0]][0][0]
                    if kiter == 1:
                        self._initial_req_rates[tuple([cts_num,cdata.index])] = req_rate
                    act_rate = sim_rates[np.where(nodes == cdata.index + 1)[0]][0]
                    if cdata.inout == "out":
                        req_tot += req_rate
                        out_tot += act_rate
                        num_ext += 1

                    # otherwise, track the inj node number for next step
                    elif cdata.inout == "in":
                        num_inj += 1
                        inj_nodes.append(cdata.index)
                        inj_req += req_rate

                if record:
                    idx = tuple([cts_num,cdata.index])
                    if idx not in self._cts_node_record_dict:
                        self._cts_node_record_dict[idx] = {"cum_req_vol":0.,"cum_act_vol":0.}
                    if cdata.inout == "out":
                        self._cts_node_record_dict[idx]["cum_act_vol"] += -act_rate * dt
                        self._cts_node_record_dict[idx]["act_rate"] = -act_rate
                    else:
                        self._cts_node_record_dict[idx]["cum_act_vol"] = 0.0
                        self._cts_node_record_dict[idx]["act_rate"] = 0.0
                    self._cts_node_record_dict[idx]["cum_req_vol"] += req_rate * dt
                    #self._cts_node_record_dict[idx]["req_rate"] = req_rate
                    self._cts_node_record_dict[idx]["req_rate"] = self._initial_req_rates[idx]

                    self._cts_node_record_dict[idx]["instance"] = cdata.instance
                    self._cts_node_record_dict[idx]["package"] = cdata.pakname
                    self._cts_node_record_dict[idx]["cts_num"] = cts_num
                    self._cts_node_record_dict[idx]["index"] = cdata.index
                    self._cts_node_record_dict[idx]["inout"] = cdata.inout

                    self._cts_current_nodes.append(idx)
                    if cdata.inout == "in":
                        inj_idx.append(idx)

            # if no extraction, return
            if req_tot == 0:
                continue

            # if non-zero extraction but no injectors, bad times
            if num_inj == 0:
                if req_tot != 0:
                    raise Exception("no injectors found and extraction is non-zero")
                else:
                    continue
            # if the requested inj rate is zero
            org_inj_req = inj_req
            if inj_req == 0:
                # if also no extraction, skip to next cts system
                if out_tot == 0:
                    continue
                # otherwise, use the simulated extraction rate
                else:
                    inj_req = out_tot

            if record:
                cts_instance.cum_req_vol += req_tot * dt
                cts_instance.cum_act_vol += out_tot * dt
                cts_instance.req_rate = req_tot
                cts_instance.act_rate = out_tot
                cts_instance.num_inj = len(inj_nodes)
                cts_instance.num_ext = num_ext

            # now set the rates for the injectors
            for cdata in cts_data:
                if cdata.index not in inj_nodes:
                    continue

                if cdata.pakname.lower() == "maw":
                    # requested flow per maw
                    addr = ["RATE", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    req_rates = self._gwf.get_value_ptr(wbaddr)
                    # the rate for this injector is a function of its requested rate
                    # vs the total inj rate
                    req_frac = req_rates[cdata.index - 1] / inj_req
                    # need to flip the sign here
                    if cts_instance._balance_flows:
                        req_rates[cdata.index - 1] = -1. * out_tot * req_frac
                    else:
                        if org_inj_req == 0:
                            req_frac = 0.0
                        else:
                            req_frac = req_rates[cdata.index - 1] / org_inj_req
                        out_tot = org_inj_req


                else:
                    # the bound vector is what we want to set - its get read into
                    # the rhs each time solve() is called
                    addr = ["BOUND", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    rates = self._gwf.get_value_ptr(wbaddr)

                    # again the node list to be careful
                    addr = ["NODELIST", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    nodes = self._gwf.get_value(wbaddr)

                    # work on the index of this occurence
                    idx = np.where(nodes == cdata.index + 1)[0][cdata.occurence]

                    # the rate for this injector is a function of its requested rate
                    # vs the total inj rate
                    #req_frac = rates[np.where(nodes == cdata.index + 1)[0]] / inj_req
                    req_frac = rates[idx] / inj_req
                    # dont need to flip the sign here bc the RHS values are already flipped
                    if cts_instance._balance_flows:
                        rates[idx] = out_tot * req_frac
                    else:
                        if org_inj_req == 0:
                            req_frac = 0
                        else:
                            req_frac = rates[idx] / org_inj_req
                        out_tot = org_inj_req
                if record:
                    idx = tuple([cts_num, cdata.index])
                    if idx in inj_idx:
                        self._cts_node_record_dict[idx]["cum_act_vol"] += float(np.abs(out_tot * req_frac * dt))
                        self._cts_node_record_dict[idx]["act_rate"] = float(np.abs(out_tot * req_frac))


    def _write_flow_summary(self, stress_period, time_step,ctime,dt):
        """write a flow summary for `stress_period`

        Parameters:
            stress_period (int): the one-based stress period to summarize

        """
        f = self._flow_node_summary

        #self._flow_node_summary.write("stress_period,time_step,ctime,dt,cts_system,package,instance,index,")
        #self._flow_node_summary.write("requested_rate,actual_rate,requested_cum_vol,actual_cum_vol\n")
        for idx in self._cts_current_nodes:
            rec = self._cts_node_record_dict[idx]
            f.write("{0},{1},{2},{3}".format(stress_period, time_step, ctime, dt))
            f.write(",{0}".format(rec["cts_num"]))
            f.write(",{0}".format(rec["package"]))
            f.write(",{0}".format(rec["instance"]))
            f.write(",{0}".format(rec["inout"]))
            f.write(",{0}".format(rec["index"]))
            f.write(",{0}".format(rec["req_rate"]))
            f.write(",{0}".format(rec["act_rate"]))
            f.write(",{0}".format(rec["cum_req_vol"]))
            f.write(",{0}\n".format(rec["cum_act_vol"]))
        if len(self._cts_current_nodes) > 0:
            f = self._flow_system_summary
            for cts_num, cts_instance in self._cts_instances.items():
                f.write("{0},{1},{2},{3}".format(stress_period, time_step, ctime, dt))
                f.write(",{0}".format(cts_num))
                f.write(",{0}".format(cts_instance.num_inj))
                f.write(",{0}".format(cts_instance.num_ext))
                f.write(",{0}".format(cts_instance.req_rate*-1.))
                f.write(",{0}".format(cts_instance.act_rate))
                f.write(",{0}".format(cts_instance.cum_req_vol*-1))
                f.write(",{0}\n".format(cts_instance.cum_act_vol))



    def solve_gwf(self):
        """solve the flow across the modflow sim times

        """
        if self._gwf is None:
            #raise Exception("gwf is None")
            self._gwf = self._initialize_gwf(self._lib_name,self._flow_dir)
        self._initial_req_rates = {}
        sim_start = datetime.now()
        print("...starting flow solution at {0}".format(sim_start.strftime(DT_FMT)))
        # reset the node tracking containers
        self._cts_node_record_dict = {}
        self._cts_current_nodes = []
        # get current sim time
        ctime = self._gwf.get_current_time()
        # get ending sim time
        etime = self._gwf.get_end_time()
        # max number of iterations
        max_iter = self._gwf.get_value(self._gwf.get_var_address("MXITER", "SLN_1"))
        # let's do it!
        num_fails = 0
        while ctime < etime:
            sol_start = datetime.now()
            # the length of this sim time
            dt = self._gwf.get_time_step()
            # prep the current time step
            self._gwf.prepare_time_step(dt)
            kiter = 0
            # prep to solve
            self._gwf.prepare_solve(1)
            # the current one-based stress period number
            stress_period = self._gwf.get_value(self._gwf.get_var_address("KPER", "TDIS"))[0]
            time_step = self._gwf.get_value(self._gwf.get_var_address("KSTP", "TDIS"))[0]
            # solve until converged
            while kiter < max_iter:
                # balance the cts flows based on current extraction rates in the RHS
                if kiter > 0:
                    self._balance_cts_flows(stress_period, dt, kiter,False)
                convg = self._gwf.solve(1)
                if convg:
                    td = (datetime.now() - sol_start).total_seconds()/60.0
                    print("flow stress period,time step {0},{1} converged with {2} iters, took {3:10.5G} mins".format(stress_period,time_step,kiter,td))
                    break

                kiter += 1

            # this time record the deets
            self._balance_cts_flows(stress_period, dt, kiter, True)
            if not convg:
                td = (datetime.now() - sol_start).total_seconds() / 60.0
                print("flow stress period,time step {0},{1} did not converge, {2} iters, took {3:10.5G} mins".format(stress_period,time_step,kiter,td))
                num_fails += 1
            try:
                self._gwf.finalize_solve(1)
            except:
                pass

            self._gwf.finalize_time_step()
            # update current sim time
            ctime = self._gwf.get_current_time()
            self._write_flow_summary(stress_period,time_step,ctime,dt)
        sim_end = datetime.now()
        td = (sim_end - sim_start).total_seconds() / 60.0
        print("\n...flow solution finished at {0}, took: {1:10.5G} mins".format(sim_end.strftime(DT_FMT),td))
        if num_fails > 0:
            print("...failed to converge {0} times".format(num_fails))
        print("\n")
        self._flow_node_summary.close()
        self._flow_system_summary.close()

    def solve_gwt(self):
        """solve transport across the modflow sim times


        """
        if self._gwt is None:
            #raise Exception("gwt is None")
            self._gwt = self._initialize_gwt(self._lib_name,self._transport_dir)

        if self._gwf is None:
            #raise Exception("gwt is None")
            self._gwf = self._initialize_gwf(self._lib_name,self._flow_dir)
        sim_start = datetime.now()
        print("...starting transport solution at {0}".format(sim_start.strftime(DT_FMT)))
        # reset the node tracking containers
        self._cts_node_record_dict = {}
        self._cts_current_nodes = []
        # get the current sim time
        ctime = self._gwt.get_current_time()
        # get the ending sim time
        etime = self._gwt.get_end_time()
        # max number of solution iterations
        max_iter = self._gwt.get_value(self._gwt.get_var_address("MXITER", "SLN_1"))
        num_fails = 0
        # let's do it!
        while ctime < etime:
            sol_start = datetime.now()
            # length of the current solve time
            dt = self._gwt.get_time_step()
            # prep the current time step
            self._gwt.prepare_time_step(dt)
            kiter = 0
            # prep to solve
            self._gwt.prepare_solve(1)
            # the one-based stress period number
            stress_period = self._gwt.get_value(self._gwt.get_var_address("KPER", "TDIS"))[0]
            time_step = self._gwt.get_value(self._gwt.get_var_address("KSTP", "TDIS"))[0]

            # solve until converged
            while kiter < max_iter:
                # apply efficiency to injectors
                self._set_current_injection_concentrations(ctime, stress_period,dt)
                convg = self._gwt.solve(1)
                if convg:
                    td = (datetime.now() - sol_start).total_seconds() / 60.0
                    print("transport stress period,time step {0},{1} converged with {2} iters, took {3:10.5G} mins".format(stress_period, time_step, kiter,td))
                    break
                kiter += 1
            # run thru one last time to record the final converged cts metrics
            self._set_current_injection_concentrations(ctime, stress_period, dt, record=True)
            if not convg:
                td = (datetime.now() - sol_start).total_seconds() / 60.0
                print("transport stress period,time step {0},{1} did not converged, {2} iters, took {3:10.5G} mins".format(
                    stress_period, time_step, kiter, td))
                num_fails += 1

            try:
                self._gwt.finalize_solve(1)
            except:
                pass
            self._gwt.finalize_time_step()

            # update the current time tracking
            ctime = self._gwt.get_current_time()
            self._write_transport_summary(stress_period,time_step,ctime,dt)
        sim_end = datetime.now()
        td = (sim_end - sim_start).total_seconds() / 60.0
        print("\n...transport solution finished at {0}, took: {1:10.5G} mins".format(sim_end.strftime(DT_FMT),td))
        if num_fails > 0:
            print("...failed to converge {0} times".format(num_fails))
        print("\n")
        self._cts_node_summary.close()
        self._cts_system_summary.close()
        if self._maw_node_summary is not None:
            self._maw_node_summary.close()

    def _write_transport_summary(self,stress_period,time_step,ctime,dt):
        f = self._cts_node_summary

        for idx in self._cts_current_nodes:
            rec = self._cts_node_record_dict[idx]
            f.write("{0},{1},{2},{3}".format(stress_period,time_step,ctime,dt))
            f.write(",{0}".format(rec["cts_num"]))
            f.write(",{0}".format(rec["package"]))
            f.write(",{0}".format(rec["instance"]))
            f.write(",{0}".format(rec["index"]))
            f.write(",{0}".format(rec["flow_rate"]))
            f.write(",{0}".format(rec["cum_vol"]))
            f.write(",{0}".format(rec["concen"]))
            f.write(",{0}".format(rec["mass_rate"]))
            f.write(",{0}".format(rec["mass"]))
            f.write(",{0}\n".format(rec["cum_mass"]))
        if len(self._cts_current_nodes) > 0:
            f = self._cts_system_summary
            for cts_num,cts_instance in self._cts_instances.items():
                eff, concen, items = cts_instance.get_eff_concen_and_entries(stress_period)
                f.write("{0},{1},{2},{3}".format(stress_period,time_step,ctime,dt))
                f.write(",{0}".format(cts_num))
                f.write(",{0}".format(cts_instance.num_inj))
                f.write(",{0}".format(cts_instance.num_ext))
                f.write(",{0}".format(cts_instance.inj_rate))
                f.write(",{0}".format(cts_instance.cum_vol))
                f.write(",{0}".format(cts_instance.mass_removed))
                f.write(",{0}".format(cts_instance.cum_mass_ext))
                f.write(",{0}".format(cts_instance.mass_treated))
                f.write(",{0}".format(cts_instance.cum_mass_treat))
                f.write(",{0}".format(cts_instance.mass_inj))
                f.write(",{0}".format(cts_instance.cum_mass_inj))
                f.write(",{0}".format(cts_instance.inj_conc))
                f.write(",{0}\n".format(eff))
        if len(self._maw_current_nodes) > 0:
            # self._maw_node_summary.write("stress_period,time_step,ctime,dt,cts_system,package,instance")
            # self._maw_node_summary.write(",maw_wellno,index,flow_rate,cum_vol,concen,mass,cum_mass\n")
            f = self._maw_node_summary
            for midx in self._maw_node_record_dict:
                rec = self._maw_node_record_dict[midx]
                f.write("{0},{1},{2},{3}".format(stress_period, time_step, ctime, dt))
                f.write(",{0}".format(rec["cts_num"]))
                f.write(",{0}".format(rec["package"]))
                f.write(",{0}".format(rec["instance"]))
                f.write(",{0}".format(rec["well_no"]))
                f.write(",{0}".format(rec["index"]))
                f.write(",{0}".format(rec["flow_rate"]))
                f.write(",{0}".format(rec["cum_vol"]))
                f.write(",{0}".format(rec["concen"]))
                f.write(",{0}".format(rec["mass_rate"]))
                f.write(",{0}".format(rec["mass"]))
                f.write(",{0}\n".format(rec["cum_mass"]))

    def _get_fmi_attr_name(self, cdata):
        """get the 'FMI-FT{} string tag for the current `CtsRecord`

        Parameters:
            cdata (CtsRecord): current cts record

        Returns:
            str: the 'FMI-FT{}' flow term tag

        """
        if cdata.pakname not in self._gwf_ft_pak_dict:
            raise Exception("CtsRecord.pakname '{0}' not in self._gwf_ft_pak_dict".format(cdata.pakname))
        return "FMI-FT{0}".format(self._gwf_ft_pak_dict[cdata.pakname])

    def _set_current_injection_concentrations(self, ctime, stress_period,dt,record=False):
        """set the injection well concentrations

        Parameters:
            ctime (float): current simulation time
            stress_period (int): one-based stress period number
            dt (float): length of the current solve
            record (bool): flag to record performance information

        """

        # get the current simulated concentration array (1-D by node number)
        conc = self._gwt.get_value_ptr(self._gwt.get_var_address("X", self._gwt_name.upper()))
        #print(conc.max())
        if record:
            self._cts_current_nodes = []
            self._maw_current_nodes = []

        # process each cts system
        for cts_num, cts_instance in self._cts_instances.items():
            cts_instance.clear_metrics()
            # get the efficiency and cts records for this stress period
            eff, concen, cts_data = cts_instance.get_eff_concen_and_entries(stress_period)
            # accumulators
            ext_mass_rate = 0.0
            inj_rate = 0.0
            ext_rate = 0.0
            inj_nodes = []
            ext_nodes = []
            inj_conc_dict = {}
            inj_idx = []
            maw_node_data = {}
            maw_inj_idx = []
            for i, cdata in enumerate(cts_data):
                # the flow rate and mass for this CtsRecord
                fr = None
                msr = None
                if cdata.pakname == "maw":
                    # mapping of gwfnode to maw wellno
                    addr = ["IMAP", self._gwf_name.upper(), cdata.instance.upper()]
                    wbaddr = self._gwf.get_var_address(*addr)
                    imap = self._gwf.get_value(wbaddr)

                    # get the FMI flow term tag
                    fmi_attr_name = self._get_fmi_attr_name(cdata)

                    # the simulated flow rates for the FMI flow term
                    addr = ["FLOW", self._gwt_name.upper(), fmi_attr_name]
                    wbaddr = self._gwt.get_var_address(*addr)
                    flow = self._gwt.get_value(wbaddr)

                    # the node numbers for this FMI flow term
                    addr = ["NODELIST", self._gwt_name.upper(), fmi_attr_name]
                    wbaddr = self._gwt.get_var_address(*addr)
                    nodes = self._gwt.get_value(wbaddr)

                    # the nodes of this maw wellno
                    maw_nodes = nodes[imap == cdata.index]

                    # the maw node flow array
                    maw_flow = np.array([flow[nodes == maw_node][0] for maw_node in maw_nodes])
                    # the  maw node conc array
                    maw_conc = np.array([conc[maw_node] for maw_node in maw_nodes])
                    # the maw mass array
                    maw_mass_rate = np.abs(maw_flow * maw_conc)
                    # get the totals
                    fr = maw_flow.sum()
                    msr = np.abs(maw_mass_rate).sum()
                    cn = maw_conc.mean()
                    maw_node_data["node"] = maw_nodes
                    maw_node_data["flow_rate"] = maw_flow
                    maw_node_data["concen"] = maw_conc
                    maw_node_data["mass_rate"] = maw_mass_rate

                else:
                    # get the FMI flow term tag
                    fmi_attr_name = self._get_fmi_attr_name(cdata)

                    # get the flow rates for this FMI flow term
                    addr = ["FLOW", self._gwt_name.upper(), fmi_attr_name]
                    wbaddr = self._gwt.get_var_address(*addr)
                    flow = self._gwt.get_value(wbaddr)

                    # get the nodes for this FMI flow term
                    addr = ["NODELIST", self._gwt_name.upper(), fmi_attr_name]
                    wbaddr = self._gwt.get_var_address(*addr)
                    nodes = self._gwt.get_value(wbaddr)

                    # one-based node compare
                    #idx = nodes == (cdata.index + 1)
                    idx = np.where(nodes==cdata.index+1)[0][cdata.occurence]

                    # the flow rate for this cts record
                    fr = flow[idx]
                    # get the node for this cts record
                    idx = nodes[idx] - 1
                    # get the conc for this cts record
                    cn = conc[idx]
                    # mass
                    msr = np.abs(fr * cn)

                if fr == 0:
                    continue

                # accumulate injectors and extractors
                if cdata.inout == "in" and fr >= 0:
                    inj_rate += fr
                    inj_nodes.append(cdata.index)
                elif cdata.inout == "out" and fr < 0:
                    ext_mass_rate += msr
                    ext_rate += fr
                    ext_nodes.append(cdata.index)
                if record:
                    idx = tuple([cts_num,cdata.index])
                    if idx not in self._cts_node_record_dict:
                        self._cts_node_record_dict[idx] = {"cum_vol":0.,"cum_mass":0.}
                    self._cts_node_record_dict[idx]["cum_vol"] += fr * dt
                    self._cts_node_record_dict[idx]["concen"] = cn
                    self._cts_node_record_dict[idx]["mass"] = msr * dt
                    self._cts_node_record_dict[idx]["mass_rate"] = msr

                    #record the ext mass and concen now - cant do inj yet...
                    if cdata.inout == "out":
                        self._cts_node_record_dict[idx]["cum_mass"] += msr * dt
                    self._cts_node_record_dict[idx]["flow_rate"] = fr
                    self._cts_node_record_dict[idx]["instance"] = cdata.instance
                    self._cts_node_record_dict[idx]["package"] = cdata.pakname
                    self._cts_node_record_dict[idx]["cts_num"] = cts_num
                    self._cts_node_record_dict[idx]["index"] = cdata.index
                    self._cts_current_nodes.append(idx)
                    if cdata.inout == "in" and fr >= 0:
                        inj_idx.append(idx)
                    if len(maw_node_data) > 0:
                        #self._maw_node_summary.write("stress_period,time_step,ctime,dt,cts_system,package,instance")
                        #self._maw_node_summary.write(",maw_wellno,index,flow_rate,cum_vol,concen,mass,cum_mass\n")
                        for node,flow,concen,mass_rate in zip(maw_node_data["node"],maw_node_data["flow_rate"],
                                                              maw_node_data["concen"],maw_node_data["mass_rate"]):
                            midx = tuple([idx,node])
                            if midx not in self._maw_node_record_dict:
                                self._maw_node_record_dict[midx] = {"cum_vol":0.0,"cum_mass":0.0}
                            self._maw_node_record_dict[midx]["well_no"] = cdata.index
                            self._maw_node_record_dict[midx]["flow_rate"] = flow
                            self._maw_node_record_dict[midx]["cum_vol"] += flow * dt
                            self._maw_node_record_dict[midx]["concen"] = concen
                            self._maw_node_record_dict[midx]["mass"] = mass_rate * dt
                            self._maw_node_record_dict[midx]["mass_rate"] = mass_rate

                            if cdata.inout == "out":
                                self._maw_node_record_dict[midx]["cum_mass"] += mass_rate * dt
                            else:
                                maw_inj_idx.append(midx)
                            self._maw_node_record_dict[midx]["instance"] = cdata.instance
                            self._maw_node_record_dict[midx]["package"] = cdata.pakname
                            self._maw_node_record_dict[midx]["cts_num"] = cts_num
                            self._maw_node_record_dict[midx]["index"] = node
                            self._maw_current_nodes.append(midx)
            # sanity checks
            if ext_mass_rate != 0 and inj_rate == 0:
                warnings.warn(
                    "injection total is zero for cts system {0}, time {1}, stress period {2}".format(cts_num, ctime,
                                                                                                     stress_period))

            # apply the efficiency factor
            if eff == 1.0:
                inj_mass_rate = 0.0
            else:
                inj_mass_rate = ext_mass_rate - (ext_mass_rate * eff)

            # calc injection concen
            inj_conc = 0.0
            if inj_rate != 0:
                inj_conc = inj_mass_rate / inj_rate
            inj_conc_dict.update({i: inj_conc for i in inj_nodes})

            if record:

                cts_instance.cum_mass_ext += ext_mass_rate * dt
                cts_instance.cum_mass_treat += ext_mass_rate * eff * dt
                cts_instance.cum_mass_inj += inj_mass_rate * dt
                cts_instance.inj_conc = inj_conc
                cts_instance.inj_rate = inj_rate
                cts_instance.mass_removed = ext_mass_rate * dt
                cts_instance.mass_treated = ext_mass_rate * eff * dt
                cts_instance.mass_inj = inj_mass_rate * dt
                cts_instance.cum_vol += inj_rate * dt
                cts_instance.num_inj = len(inj_nodes)
                cts_instance.num_ext = len(ext_nodes)



            # set the concentration for the injectors
            for i, cdata in enumerate(cts_data):
                if cdata.pakname == "maw":
                    # get the CONCRATE pointer - this is where the MWT input
                    # concentration is stored
                    # assuming there is only one MAW and only one MWT
                    # todo: support multiple maw-mwt pairs
                    if cdata.index in inj_nodes:
                        # addr = ["CONCBUDSSM", self._gwt_name.upper(), self._mwt_name.upper()]
                        # wbaddr = self._gwt.get_var_address(*addr)
                        # t1 = self._gwt.get_value_ptr(wbaddr)
                        #
                        # addr = ["CONCFEAT", self._gwt_name.upper(), self._mwt_name.upper()]
                        # wbaddr = self._gwt.get_var_address(*addr)
                        # t2 = self._gwt.get_value_ptr(wbaddr)

                        addr = ["CONCRATE", self._gwt_name.upper(), self._mwt_name.upper()]
                        wbaddr = self._gwt.get_var_address(*addr)
                        cts_conc = self._gwt.get_value_ptr(wbaddr)
                        #print(inj_conc)
                        cts_conc[cdata.index - 1] = inj_conc


                else:
                    # get the FMI flow term tag
                    fmi_attr_name = self._get_fmi_attr_name(cdata)
                    # get the nodes for this FMI flow term
                    addr = ["NODELIST", self._gwt_name.upper(), fmi_attr_name]
                    wbaddr = self._gwt.get_var_address(*addr)
                    nodes = self._gwt.get_value(wbaddr)

                    # get the AUXVAR pointer for this FMI flow term - this is where concentration is stored
                    addr = ["AUXVAR", self._gwt_name.upper(), fmi_attr_name]
                    wbaddr = self._gwt.get_var_address(*addr)
                    cts_conc = self._gwt.get_value_ptr(wbaddr)

                    # for each injector, set the concentration
                    for inj_node, inj_conc in inj_conc_dict.items():
                        idx = np.where(nodes == inj_node + 1)[0]
                        cts_conc[idx] = inj_conc

            #now record inj metrics
            if record:
                for idx in inj_idx:
                    self._cts_node_record_dict[idx]["cum_mass"] += inj_conc * self._cts_node_record_dict[idx]["flow_rate"] * dt
                    self._cts_node_record_dict[idx]["concen"] = inj_conc
                    self._cts_node_record_dict[idx]["mass"] = inj_conc * self._cts_node_record_dict[idx]["flow_rate"] * dt
                    self._cts_node_record_dict[idx]["mass_rate"] = inj_conc * self._cts_node_record_dict[idx]["flow_rate"]
                for midx in maw_inj_idx:
                    self._maw_node_record_dict[midx]["cum_mass"] += inj_conc * self._maw_node_record_dict[midx]["flow_rate"] * dt
                    self._maw_node_record_dict[midx]["concen"] = inj_conc
                    self._maw_node_record_dict[midx]["mass"] = inj_conc * self._maw_node_record_dict[midx]["flow_rate"] * dt
                    self._maw_node_record_dict[midx]["mass"] = inj_conc * self._maw_node_record_dict[midx]["flow_rate"]

        return

    def _read_cts_input_file(self):
        """read the cts input file

        """
        # used to detect location-pak duplicates
        current_period = -1
        current_entries = []

        with open(self.cts_filename, 'r') as f:
            count = 0
            while True:
                line = f.readline()
                count += 1
                # eof
                if line == "":
                    break

                # skip empty lines or comment lines
                if len(line.strip()) == 0 or line.strip()[0] == "#":
                    continue

                # read the options block
                if line.lower().strip().startswith("begin options"):
                    while True:
                        line2 = f.readline()
                        count += 1

                        if line2 == "":
                            raise EOFError("EOF while reading options")
                        elif len(line.strip()) == 0 or line.strip()[0] == "#":
                            continue
                        elif line2.lower().strip().startswith("begin"):
                            raise Exception("a new begin block found while parsing options")
                        elif line2.lower().strip().startswith("end options"):
                            break

                        # todo: parse options


                # parse a new cts system period block
                elif line.lower().strip().startswith("begin period"):
                    raw = line.lower().strip().split()
                    try:
                        period = int(raw[2])
                    except:
                        raise Exception("error parsing period number on line {0} ('{1}')".format(count, line))
                    try:
                        cts_system_num = int(raw[4])
                    except:
                        raise Exception("error parsing cts system number on line {0} ('{1}')".format(count, line))
                    # default efficency = 1.0
                    eff = 1.0
                    if "efficiency" in line.lower():
                        try:
                            #eff = float(raw[-1])
                            eff = float(line.lower().split("efficiency")[1].split()[0])
                        except:
                            raise Exception("error parsing efficiency on line {0} ('{1}')".format(count, line))
                    if eff > 1.0 or eff < 0.0:
                        raise Exception("invalid efficiency {0} on line {1} ('{2}')".format(eff,count,line))

                    concen = 0.0
                    if "concentration" in line.lower():
                        raise NotImplementedError("concentration not implemented")
                    if concen < 0.0:
                        raise Exception("invalid concentration {0} on line {1} ('{2}')".format(concen,count,line))

                    #for tracking multiple occurences in the same cell
                    if period != current_period:
                        current_period = period
                        current_entries = []

                    cts_system_entries = []
                    num_inj, num_ext = 0,0
                    while True:
                        line2 = f.readline()
                        count += 1
                        if line2 == "":
                            raise EOFError("EOF while reading cts system block '{0}'".format(line))
                        elif len(line.strip()) == 0 or line.strip()[0] == "#":
                            continue
                        elif line2.lower().strip().startswith("begin"):
                            raise Exception("a new begin block found while parsing cts system block '{0}'".format(line))
                        elif line2.lower().strip().startswith("end period"):
                            break

                        raw = line2.lower().strip().split()
                        inout = raw[2].lower()
                        if inout not in ["in", "out"]:
                            raise Exception("third entry on cts line {0} ('{1}') should be 'in' or 'out', not '{2}'". \
                                            format(count, line2, inout))
                        if inout == "in":
                            num_inj += 1
                        else:
                            num_ext += 1

                        bc_name = raw[1].lower()
                        package_type = raw[0].lower()
                        # if package_type not in self.cts_bc_texts:
                        #     self.cts_bc_texts.append(package_type)

                        # if maw, read wellno
                        if package_type == "maw":
                            wellno = int(raw[3])
                            cts_system_entries.append([wellno, inout, package_type, bc_name])

                        # else if structured, read l-r-c
                        elif self.is_structured:
                            if len(raw) < 6:
                                raise Exception("wrong number of cts entries for structured grid on line {0}: '{1}' ". \
                                                format(count, line2))
                            kij = []
                            for i in range(3):
                                try:
                                    kij.append(int(raw[i + 3]) - 1)
                                except:
                                    raise Exception("error casting k-i-j info on line {0}: '{1}'".format(count, line2))
                            # todo: check that the kij info is inbounds

                            # convert to node number
                            cts_system_entries.append(
                                [self._structured_mg.get_node([kij])[0], inout, package_type, bc_name])

                        else:
                            raise NotImplementedError("only structured grids currently supported")
                    if len(cts_system_entries) != 0:
                        #warnings.warn("no entries found for cts system {0}".format(cts_system_num))
                        if num_inj == 0:
                            raise Exception("no injectors found for period {0} cts system {1}".format(period,cts_system_num))
                        elif num_ext == 0:
                            raise Exception(
                                "no extractors found for period {0} cts system {1}".format(period, cts_system_num))

                    if len(cts_system_entries) == 0:
                        print("...cts system {0} is explicitly 'off' in period {1}".format(cts_system_num,period))
                    # if this is a new cts system, instaniate
                    if cts_system_num not in self._cts_instances:
                        self._cts_instances[cts_system_num] = CtsSystem(cts_system_num, period, eff, concen,
                                                                        cts_system_entries)
                    # otherwise, add entries
                    else:
                        self._cts_instances[cts_system_num].add_period_entries(period, eff, concen,
                                                                               cts_system_entries)
                    for cts_rec in self._cts_instances[cts_system_num]._entries[period]:
                        t = cts_rec.get_tuple_wo_instance()
                        if t in current_entries:
                            raise Exception("duplicate locations across cts instances detected for period"+
                                            " {0} - this currently not supported".format(period))
                        current_entries.append(t)

                else:
                    raise Exception("unrecognized cts file input on line {0}: '{1}'".format(count, line))
        if len(self._cts_instances) == 0:
            raise Exception("no cts systems found in cts file")


def main():
    """command line driver

    Note:
        command line usage requires a config python source file.  This file must have all the
        args required to instantiate Mf6Cts as well as the flow model output files that
        the FMI package needs to run (these should be in a python list structure)

    """
    args = sys.argv
    if len(args) != 2:
        print(args)
        raise Exception("mf6cts commndline usage: 'python mf6cts.py <config_file.py>")
    config_file = args[1]
    if not os.path.exists(config_file):
        raise Exception("mf6cts config_file '{0}' not found".format(config_file))
    config_module = __import__(config_file.replace(".py",""))

    mf6 = Mf6Cts(config_module.cts_filename, config_module.lib_name, config_module.transport_dir, config_module.flow_dir,
                 config_module.is_structured)
    if hasattr(config_module,"balance_flows"):
        for cts_num,cts_instance in mf6._cts_instances.items():
            cts_instance._balance_flows = config_module.balance_flows
    solve_gwf = True
    if hasattr(config_module,"solve_gwf"):
        solve_gwf = bool(config_module.solve_gwf)
    if solve_gwf:
        mf6.solve_gwf()

    transfer_flow_output_files = True
    if hasattr(config_module,"transfer_flow_output_files"):
        transfer_flow_output_files = bool(config_module.transfer_flow_output_files)
    if transfer_flow_output_files:
        for flow_output_file in config_module.flow_output_files:
            shutil.copy2(os.path.join(config_module.flow_dir, flow_output_file),
                         os.path.join(config_module.transport_dir, flow_output_file))

    solve_gwt = True
    if hasattr(config_module, "solve_gwt"):
        solve_gwt = bool(config_module.solve_gwt)
    if solve_gwt:
        mf6.solve_gwt()

    mf6.finalize()


if __name__ == "__main__":
    main()