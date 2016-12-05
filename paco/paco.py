#!/usr/bin/env python
import logging
import optparse
import os.path
import pacoRun
import pacoUtils as pu
import yaml

stepsToRun = [pacoRun.runStage1, pacoRun.runStage2, pacoRun.runStage3]

def main():
    p = optparse.OptionParser(
            usage  ="usage: %prog [options] CONFIG [DIRECTORY]",
            description =   "Execute run with config file."
                            "Optionally, can specify a directory"
                            "(if none given, uses same as config file.)",
            version="%prog 2.0"
            )
    p.add_option('--prefix', '-p',
            help="override prefix with which all intermediate and output files with this prefix.")
    p.add_option('--verbose', '-v', action="store_true", default=False)
    CLIoptions, arguments = p.parse_args()
    if len(arguments) > 2 or len(arguments) < 1:
        p.error("incorrect number of arguments, must specify the config file (and optionally a directory).")

    # set up logging
    logLevel = logging.WARNING
    if CLIoptions.verbose:
        logLevel = logging.DEBUG

    logging.basicConfig(
            level=logLevel,
            format='%(asctime)s. %(levelname)s in %(funcName)s: %(message)s')

    # get first arg: configPath
    configPath = arguments[0]

    # set up paths
    configDir, configBasename = os.path.split(configPath)
    configDir = os.path.abspath(configDir)

    # get second arg: dataDir
    if len(arguments) > 1:
        dataDir = arguments[1]
    else:
        dataDir = configDir
    dataDir = os.path.abspath(dataDir)

    # if not os.path.isdir(dataDir):
    #     msg = "directory not found: %s" % dataDir
    #     logging.error(msg)
    #     return -1

    pu.ensurePath(dataDir)

    # now, read config
    # if the file does not end with yaml, warn
    fileExtension = os.path.splitext(configBasename)[1]
    if fileExtension != ".yaml":
        logging.warn("config file does not end with .yaml: %s", configPath)

    msg = "Executing config file %s" % configPath
    print(msg)
    logging.info(msg)

    with open(configPath, 'r') as f:
        params = yaml.load(f)

    # print(params)

    # add/override options in params (prefix, dataDir, etc.)
    params['Config']['verbose'] = CLIoptions.verbose
    if CLIoptions.prefix is None:
        if 'Prefix_Default' in params['Config']:
            params['Config']['Prefix'] = params['Config']['Prefix_Default']
        else:
            logging.error("No prefix provided, and no Prefix_Default in config file.")
            return -1
    else:
        if 'Prefix_Default' in params['Config']:
            msg = "Prefix provided (%s) used instead of Prefix_Default (%s)" % (
                CLIoptions.prefix, params['Config']['Prefix_Default']
            )
            print(msg)
            logging.warn(msg)
        params['Config']['Prefix'] = CLIoptions.prefix

    # record directories
    pacoRootDir = os.path.abspath(os.path.dirname(pacoRun.__file__))
    params['Config']['DirData'] = dataDir
    params['Config']['DirCurrent'] = os.getcwd()
    params['Config']['DirConfig'] = configDir
    params['Config']['DirPacoRoot'] = pacoRootDir
    params['Config']['DirPacoBin'] = os.path.abspath(os.path.join(pacoRootDir, '..', 'src'))

    # what about stokes.edp?    # deal with in runStage1

    # dump actually used params to new location
    pu.ensurePath(dataDir)
    if not configBasename.startswith(params['Config']['Prefix']):
        configBasename_new = params['Config']['Prefix'] + "s0_" + configBasename
    else:
        configBasename_new = configBasename
    with open(os.path.join(dataDir, configBasename_new), "w") as f:
        yaml.dump(params, f)

    runStages(params)
    logging.shutdown()


def runStages(params):
    # Now run the stuff:
    for func in stepsToRun:
        fname = func.func_name
        if fname in params:
            msg = "Executing %s." % fname
            print(msg)
            logging.info(msg)
            func(params)


if __name__ == '__main__':
    main()
