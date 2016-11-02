#!/usr/bin/env python
import logging
import optparse
import os.path
import pacoHelpers


def main():
    p = optparse.OptionParser(
            usage  ="usage: %prog [options] sourcePath",
            version="%prog 1.0"
            )
    p.add_option('--prefix', '-p', default="pre_",
            help="prefix all intermediate and output files with this prefix")
    p.add_option('--destPath', '-d')
    p.add_option('--verbose', '-v', action="store_true", default=False)
    options, arguments = p.parse_args()
    if len(arguments) != 1:
        p.error("incorrect number of arguments")

    sourcePath = arguments[0]
    sourceDir = os.path.dirname(sourcePath)

    logLevel = logging.WARNING
    if options.verbose:
        logLevel = logging.DEBUG

    logging.basicConfig(
            level=logLevel,
            format='%(asctime)s. %(levelname)s in %(funcName)s: %(message)s')

    if options.destPath is not None:
        pass

    logging.warn("Running S1 on %s with prefix %s."
                    % (sourcePath, options.prefix))
    pacoHelpers.runS1(sourcePath=sourcePath, prefix=options.prefix)
    logging.warn("Running S2 on %s with prefix %s."
                    % (sourceDir, options.prefix))
    pacoHelpers.runS2(sourcePath=sourceDir, prefix=options.prefix)

    logging.shutdown()

if __name__ == '__main__':
    main()
