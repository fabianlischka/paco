#!/usr/bin/env python
import logging
import optparse
import os.path
import pacoHelpers


def main():
    p = optparse.OptionParser(
            usage  ="usage: %prog [options] sourceFullPath",
            description =   "Call with a FreeFem++ *.edp file as argument. "
                            "It will then run FreeFem++, FVM, etc. on it.",
            version="%prog 1.0"
            )
    p.add_option('--prefix', '-p', default="pre_",
            help="prefix all intermediate and output files with this prefix")
    p.add_option('--destDir', '-d',
            help="store all intermediate and output files in this directory")
    p.add_option('--verbose', '-v', action="store_true", default=False)
    options, arguments = p.parse_args()
    if len(arguments) != 1:
        p.error("incorrect number of arguments, must specify the sourceFullPath.")

    sourceFullPath = arguments[0]

    runStages(sourceFullPath, options.prefix, options.destDir, options.verbose)


def runStages(sourceFullPath, prefix, destDir, verbose):
    # set up logging
    logLevel = logging.WARNING
    if verbose:
        logLevel = logging.DEBUG

    logging.basicConfig(
            level=logLevel,
            format='%(asctime)s. %(levelname)s in %(funcName)s: %(message)s')

    # set up paths
    sourceDir, sourceBasename = os.path.split(sourceFullPath)
    if destDir is None:
        destDir = sourceDir

    # S1
    msg = "Run stage 1 on %s with prefix %s into %s" % (sourceFullPath, prefix, destDir)
    print(msg)
    logging.info(msg)
    pacoHelpers.runS1(sourceFullPath=sourceFullPath, prefix=prefix, destDir=destDir)
    # S2
    msg = "Run stage 2 with prefix %s in %s." % (prefix, destDir)
    print(msg)
    logging.info(msg)
    pacoHelpers.runS2(sourceDir=destDir, prefix=prefix)

    logging.shutdown()





if __name__ == '__main__':
    main()
