import collections
import ROOT

class keydefaultdict(collections.defaultdict):
    """
    http://stackoverflow.com/a/2912455
    """
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret
tfiles = keydefaultdict(ROOT.TFile.Open)

def run1or2(*args):
    if True:#"ggH" in args:
        return 2
    else:
        return 1


def gettemplate(*args):
    run = run1or2(*args)
    print "run", run, args
    argsbkp = args
    args = list(args)
    productionmode = None
    hypothesis = None
    flavor = None

    for a in "ggH", "ggZZ", "qqZZ", "ZX", "data":
        if a in args:
            productionmode = a
            args.remove(a)
            break
    if productionmode is None:
        raise ValueError("No productionmode in {}".format(argsbkp))

    if productionmode == "ggH":
        for a in "0+", "0-", "int":
            if a in args:
                hypothesis = a
                args.remove(a)
                break
        if hypothesis is None:
            raise ValueError("No hypothesis in {}".format(argsbkp))

    for a in "2e2mu", "4e", "4mu":
        if a in args:
            flavor = a
            args.remove(a)
            break
    if flavor is None:
        raise ValueError("No flavor in {}".format(argsbkp))

    if args:
        raise ValueError("Extra arguments in {}".format(argsbkp))

    isbkg = (productionmode != "ggH")
    if productionmode == "data":
        filename = "templates_Heshy/data.root"
    else:
        if run == 2:
            filename = "templates_Heshy/{}{}_fa3Adap_new.root".format(flavor, "_bkg" if isbkg else "")
        elif run == 1:
            filename = "templates_run1/{}_fa3Adap_new{}.root".format(flavor, "_bkg" if isbkg else "")
    if productionmode == "ggH":
        templatename = {
                        "0+": "template0PlusAdapSmoothMirror",
                        "0-": "template0MinusAdapSmoothMirror",
                        "int": "templateg1g4AdapSmooth",
                       }[hypothesis]
    elif productionmode == "data":
        templatename = "candTree"
    else:
        if run == 2:
            templatename = "template{}AdapSmoothMirror".format(productionmode)
        elif run == 1:
            templatename = "template_{}".format(productionmode)

    f = tfiles[filename]
    if not f:
        raise OSError("No file {}".format(f))
    h = getattr(f, templatename)
    if not h:
        raise OSError("No template {} in {}".format(templatename, filename))
    return h
