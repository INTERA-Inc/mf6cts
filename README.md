# mf6cts: Contaminant treatment system modeling via the MODFLOW 6 API

The contaminant treatment system (CTS) package is contained in the `\mf6cts\mf6cts.py:Mf6Cts` class.  Several examples usages can be found in the `autotest\cts_mf6_test.py` script - these tests are "self-contained", in that the precompiled MODFLOW 6 binary library files, as well as the required python dependencies are all present in this repo, so one should be able to "just run" these tests ("should")...

## testing


![CI workflow](https://github.com/INTERA-Inc/mf6cts/actions/workflows/ci.yml/badge.svg)


## basic user information

The `Mf6Cts` class is designed to work with separate flow-then-transport
modeling scheme. As such, the flow model and transport models should be
in separate directories. Additionally, the `FMI` process is expected to
be used, and, if the `MAW` package is used in the flow model, then the
`MWT` package is expected in the transport model. Users are expected to
make the flow-model binary output files needed by the FMI available.

The CTS input file follows the standard MODFLOW 6 input format style.
Through the use of period blocks, a CTS system can be (re-)defined for
each stress period in the model. Each period block must identify the CTS
system that is being defined. In this way, numerous different CTS
systems can be (re-)defined at the stress period level. An optional
efficiency can be specified; efficiency ranges from 0.0 (no
dissolved-phase mass is removed) to 1.0 (all dissolved-phase mass is
removed). If not specified, efficiency is assumed 1.0---the CTS instance
removes all mass.

Within each period block, information related to the components of the
CTS must be defined. These include the component package type (e.g.
"wel", "maw", etc), the component instance of the package type (as
defined in the flow model nam file, e.g. "wel_1", "maw_cts", etc),
whether the component is an extractor ("out") or an injector ("in"), and
the indexing information for the component. For non-MAW type components,
if the model is a structured grid, then the indexing information is the
layer, row and column location of the component; if the model is
unstructured, then the indexing information should be node-number of the
component. For MAW-type components, the indexing information is the
"wellno" value in the MAW package.

Below is an example CTS input file:

`begin period 2 cts 1 efficiency 0.361 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 out 3 21 21 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 out 3 21 23 `   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`maw maw_0 out 2 `   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 in 1 1 1 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`maw maw_0 in 1 `   
`end period 2 cts 1 `  
  
`begin period 2 cts 2 efficiency 0.909 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 out 3 20 3 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 in 1 33 1 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 in 1 33 33 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 in 1 1 33 `  
`end period 2 cts 2 `  
  
`begin period 5 cts 1 efficiency 0.925 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`maw maw_0 out 2 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`maw maw_0 in 1 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`end period 3 cts 1 `  
  
`begin period 7 cts 2 efficiency 0.222 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 out 3 20 3 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 in 1 33 1 `  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 in 1 33 33 `   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`wel wel_0 in 1 1 33 `  
`end period 3 cts 2 `  
  

In this example file, we see that for stress period 2, CTS instance "1"
has three extraction wells and one injection well, while CTS instance
"2" has one extraction well and three injection wells. Using the MODFLOW
6 convention, CTS instance "1" continues from stress period 2 to stress
period 5 as it is defined in stress period 2. Similarly, CTS instance 2
continues from stress period 2 to stress period 7 as it was defined in
stress period 2. At stress period 5, CTS instance "1" is reconfigured to
be a single injector and single extractor system of just MAW-type wells.
At stress period 7, CTS instance "2" changes its efficiency value but
contains the same injection and extraction wells.

The `Mf6Cts` class can also be driven from the command line:\
`python mf6cts.py <config_file.py> `\
where `<config_file.py>` is the name of a configuration file that is a
simple python source file listing the required arguments:

-   `cts_filename`: the name of the cts file

-   `lib_name`: the name of the MF6 shared library file

-   `transport_dir`: the directory holding the transport model files

-   `flow_dir`: the directory holding the flow model files

-   `is_structured`: a boolean flag indicating if the model is a
    structured grid (a value of `1` indicates a structured grid)

-   `flow_output_files`: a python list of the flow model output binary
    files that the `FMI` package is expecting

An example configuration file:

`cts_filename=’model.cts’`   
`lib_name=’libmf6.so’ `    
`transport_dir=’fivespot_maw_t_api’ `  
`flow_dir=’fivespot_maw_api’ `   
`is_structured=True `  
`flow_output_files=[’gwf.hds’,’gwf.bud’,’gwf.maw.bud’] `  
  

The MODFLOW 6 CTS package write several comma-separated-value (CSV)
files summarizing the performance of the CTS system for the both the
flow and transport models and from both the CTS node and CTS system
perspective. The flow model summary CSV files are written in the flow
model directory (e.g. `flow_dir`) and are named
"gwf_cts_flow_node_summary.csv" and "gwf_cts_flow_system_summary.csv";
these files summarize the flow-model extraction and injection aspects of
CTS instance performance and include information about the requested and
actual rates and cumulative volumes extracted and/or injected for each
CTS system and its nodes across all flow solution stress periods and
time steps.

The transport-model summary CSV files are written in the transport model
directory (e.g. `transport_dir`) and are named
"gwt_cts_node_summary.csv" and "gwt_cts_system_summary.csv"; these files
summarize the transport-model extraction and injection aspects of CTS
instance performance and include information about the extracted and
injected dissolved-phase mass as both a rate and cumulative mass for
each stress period and time step, as well as cumulative mass across the
entire solution period. These files also contain information about the
blended injected concentration if an efficiency less than 1.0 is used.

If MAW-type boundary conditions are included in CTS instance, then a
MAW-specific summary file is written in the transport model directory
named "gwt_maw_node_summary.csv". This file summarizes the
dissolved-phase transport aspects of individual nodes that comprise each
MAW boundary condition, including flow rate, concentration, cumulative
volume and cumulative mass for each node in each MAW boundary condition
that is included in a CTS instance.
