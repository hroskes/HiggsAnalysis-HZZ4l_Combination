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
    if not ("ZX" in args and "4mu" in args):
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
    systematic = None

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

        for a in "ScaleUp", "ScaleDown", "ResUp", "ResDown", "ScaleResUp", "ScaleResDown":
            if a in args:
                systematic = a
                args.remove(a)
                break

    if productionmode == "ZX":
        for a in "ZXUp", "ZXDown":
            if a in args:
                systematic = a
                args.remove(a)
                break

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
    appendname = ""
    if isbkg:
        appendname += "_bkg"
    if systematic is not None and not (run==1 and productionmode == "ZX"):
        appendname += "_"+systematic

    if productionmode == "data":
        filename = "templates_run2/data.root"
    else:
        if run == 2:
            filename = "templates_run2/{}{}_fa3Adap_new.root".format(flavor, appendname)
        elif run == 1:
            filename = "templates_run1/{}_fa3Adap_new{}.root".format(flavor, appendname)
    if productionmode == "ggH":
        templatename = {
                        "0+": "template0PlusAdapSmoothMirror",
                        "0-": "template0MinusAdapSmoothMirror",
                        "int": "templateIntAdapSmoothMirror",
                       }[hypothesis]
    elif productionmode == "data":
        templatename = "candTree"
    else:
        if run == 2:
            templatename = "template{}AdapSmoothMirror".format(productionmode)
        elif run == 1:
            if productionmode == "ZX" and systematic == "ZXUp":
                templatename = "template_qqZZ"
            elif productionmode == "ZX" and systematic == "ZXDown":
                templatename = "T_mirror"
            else:
                templatename = "template_{}".format(productionmode)

    f = tfiles[filename]
    if not f:
        raise OSError("No file {}".format(f))
    try:
        h = getattr(f, templatename)
        if not h:
            raise AttributeError
    except AttributeError:
        f.ls()
        raise OSError("No template {} in {}".format(templatename, filename))
    return h
