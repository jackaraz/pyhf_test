#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json, os, argparse, sys
from utils import HFsignal, HFbackground, pyhf_sig95Wrapper, pyhf_wrapper

parser = argparse.ArgumentParser(
    description="Run pyhf on Ma5 output."
)


modelgroup = parser.add_argument_group("Input handling")

modelgroup.add_argument("-xsec", "--xsection", default=0.000781, dest="XSEC", type=float,
                        help="Cross section for the process. Default is 0.000781 pb")
modelgroup.add_argument("--sig-path", default="./", dest="SIGPATH", type=str,
                        help="Histfactory files path for the signal. Default ./")
modelgroup.add_argument(
    "--bkg-path",
    default="./",
    dest="BKGPATH", type=str,
    help="Histfactory files path for the background.Default ./ "
)
modelgroup.add_argument("--pyhf-path", default="./tools/pyhf/src", dest="PYHFPATH", type=str,
                        help="pyhf path. Default (assuming program runs under ma5 main folder) "
                             "./tools/pyhf/src")

args = parser.parse_args()

import sys
sys.path.insert(0, args.PYHFPATH)

signal_file_names = [
    "RegionA_sig.json","RegionB_sig.json","RegionC_sig.json"
]
background_file_names = [
    "atlas_susy_2018_31_SRA.json","atlas_susy_2018_31_SRB.json","atlas_susy_2018_31_SRC.json"
]

likelihood_profile_regs = ["RegionA","RegionB","RegionC"]

for fl in signal_file_names:
    if fl not in os.listdir(args.SIGPATH):
        raise Exception(f"Cant find signal files in {args.SIGPATH}")
for fl in background_file_names:
    if fl not in os.listdir(args.BKGPATH):
        raise Exception(f"Cant find bkg files in {args.BKGPATH}")

with open(os.path.join(args.SIGPATH, "regiondata.json"), "r") as f:
    regiondata = json.load(f)
regiondata['pyhf'] = {}

minsig95 = 1e99
bestreg  = []
for fsig, fbkg, profileID in zip(signal_file_names, background_file_names, likelihood_profile_regs):

    with open(os.path.join(args.SIGPATH,fsig),"r") as f:
        signal_patch = json.load(f)
    with open(os.path.join(args.BKGPATH,fbkg),"r") as f:
        bkg_profile = json.load(f)

    HF_signal     = HFsignal(signal_patch, args.XSEC)
    HF_background = HFbackground(bkg_profile)

    regiondata['pyhf'][profileID]           = {}
    regiondata['pyhf'][profileID]["s95exp"] = -1
    regiondata['pyhf'][profileID]["s95obs"] = -1

    regiondata = pyhf_sig95Wrapper(HF_signal, HF_background, regiondata, profileID, "exp")
    regiondata = pyhf_sig95Wrapper(HF_signal, HF_background, regiondata, profileID, "obs")

    CLs = pyhf_wrapper(HF_background.hf, HF_signal())
    CLs_out = CLs['CLs_obs']
    regiondata['pyhf'][profileID]['full_CLs_output'] = CLs
    if CLs_out >= 0.:
        regiondata['pyhf'][profileID]['CLs']  = CLs_out
    s95 = float(regiondata['pyhf'][profileID]['s95exp'])
    if 0. < s95 < minsig95:
        regiondata['pyhf'][profileID]["best"] = 1
        for mybr in bestreg:
            regiondata['pyhf'][mybr]["best"]=0
        bestreg = [profileID]
        minsig95 = s95
    else:
        regiondata['pyhf'][profileID]["best"]=0

output = open("output.json","w")
output.write(json.dumps(regiondata, indent=4))
output.close()
