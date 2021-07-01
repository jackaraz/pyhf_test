import copy

class HFsignal:
    def __init__(self,HFdat, xsec):
        self.hf = HFdat
        self.xsec = xsec

    def __call__(self, xsec=None):
        HF = copy.deepcopy(self.hf)

        if xsec is None:
            xsec = self.xsec

        for iHF, inpt in enumerate(HF):
            for idat , dat in enumerate(inpt["value"]["data"]):
                HF[iHF]["value"]["data"][idat] *= xsec

        return HF



class HFbackground:
    def __init__(self,HFdat):
        self.hf = HFdat

    def impose_expected(self):
        """
            To switch observed data with total expected data per SR bin.
        """
        total_expected = {}
        HF             = copy.deepcopy(self.hf)
        for i in range(len(HF.get('observations',[]))):
            total_expected[HF['observations'][i]['name']] = [0.0]*len(HF['observations'][i]['data'])

        for iSR in range(len(HF['channels'])):
            for sample in range(len(HF['channels'][iSR]['samples'])):
                for SRbin in range(len(HF['channels'][iSR]['samples'][sample]['data'])):
                    total_expected[HF['channels'][iSR]['name']][SRbin] += \
                        HF['channels'][iSR]['samples'][sample]['data'][SRbin]

        # replace the observed bkg with total expected bkg
        for key, item in total_expected.items():
            for iobs in range(len(HF['observations'])):
                if key == HF['observations'][iobs]['name']:
                    HF['observations'][iobs]['data'] = [round(x,5) for x in item]

        return HF

    def get_expected(self):
        return self.impose_expected().get('observations',[])

    def get_observed(self):
        return self.hf.get('observations',[])



def pyhf_wrapper(background, signal):
    import pyhf
    print("pyhf has been imported from "+" ".join(pyhf.__path__))
    from pyhf.optimize import mixins

    pyhf.set_backend('numpy', precision="64b")

    try:
        workspace = pyhf.Workspace(background)
        model     = workspace.model(patches=[signal],
                                    modifier_settings={'normsys': {'interpcode': 'code4'},
                                                       'histosys': {'interpcode': 'code4p'}})
    except (pyhf.exceptions.InvalidSpecification, KeyError) as err:
        print("Invalid JSON file!! "+str(err))
        return {'CLs_obs':-1. , 'CLs_exp' : -1.}
    except Exception as err:
        print(signal)
        print("Unknown error, check pyhf_wrapper_py3 "+ str(err))
        return {'CLs_obs':-1. , 'CLs_exp' : -1.}

    def get_CLs(**kwargs):
        try:
            CLs_obs, CLs_exp = pyhf.infer.hypotest(kwargs.get('mu',1.),
                                                   workspace.data(model),
                                                   model,
                                                   test_stat="qtilde",
                                                   par_bounds=kwargs.get('bounds',
                                                                         model.config.suggested_bounds()),
                                                   return_expected=True)

        except (AssertionError, pyhf.exceptions.FailedMinimization) as err:
            print(str(err))
            # dont use false here 1.-CLs = 0 can be interpreted as false
            return 'update bounds'

        return {'CLs_obs':1.-CLs_obs , 'CLs_exp' : 1.- CLs_exp}

    #pyhf can raise an error if the poi_test bounds are too stringent
    #they need to be updated dynamically.
    update_bounds = model.config.suggested_bounds()
    iteration_limit = 0
    while True:
        CLs = get_CLs(bounds=update_bounds)
        if CLs == 'update bounds':
            update_bounds[model.config.poi_index] = (0,2*update_bounds[model.config.poi_index][1])
            iteration_limit += 1
        elif isinstance(CLs, dict):
            break
        else:
            iteration_limit += 1
        # hard limit on iteration required if it exceeds this value it means
        # Nsig >>>>> Nobs
        if iteration_limit>=3:
            return {'CLs_obs':1. , 'CLs_exp' : 1.}

    return CLs



def pyhf_sig95Wrapper(HF_signal, background, regiondata, likelihood_profile, tag):
    background = background.impose_expected() if tag == "exp" else background.hf

    def sig95(xsection):
        signal = HF_signal(xsec=xsection)
        return pyhf_wrapper(background, signal)['CLs_'+tag]-0.95

    low, hig = 1., 1.;
    while pyhf_wrapper(background, HF_signal(xsec=low))['CLs_'+tag] > 0.95:
        low *= 0.1
    while pyhf_wrapper(background, HF_signal(xsec=hig))['CLs_'+tag] < 0.95:
        hig *= 10.
    try:
        import scipy
        s95 = scipy.optimize.brentq(sig95,low,hig,xtol=low/100.)
    except:
        print("problem with scipy")
        s95=-1
    regiondata['pyhf'][likelihood_profile]["s95"+tag] = "{:.7f}".format(s95)
    print(f"{likelihood_profile}: {regiondata['pyhf'][likelihood_profile]}")

    return regiondata