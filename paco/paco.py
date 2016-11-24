#!/usr/bin/env python
import logging
import optparse
import os.path
import pacoRun
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
    if options.verbose:
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

    if not os.path.isdir(dataDir):
        msg = "directory not found: %s" % dataDir
        logging.error(msg)
        return -1

    dataDir = os.path.abspath(dataDir)

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

    # add/override prefix, dataDir in params

    # what about stokes.edp?

    # dump actually used params to new location


    runStages(params)


def runStages(params):
    # Now run the stuff:
    for func in stepsToRun:
        fname = func.func_name
        if fname in params:
            msg = "Executing %s." % fname
            func(params)


    # S1
    msg = "Run stage 1 on '%s' with prefix '%s' into dir '%s'." % (sourceFullPath, prefix, destDir)
    print(msg)
    logging.info(msg)
    pacoRun.runS1(sourceFullPath=sourceFullPath, prefix=prefix, destDir=destDir)
    # S2
    msg = "Run stage 2 with prefix '%s' in dir '%s'." % (prefix, destDir)
    print(msg)
    logging.info(msg)
    pacoRun.runS2(sourceDir=destDir, prefix=prefix)

    logging.shutdown()





if __name__ == '__main__':
    main()
