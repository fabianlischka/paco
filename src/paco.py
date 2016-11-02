#!/usr/bin/env python
import optparse
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
    print("Running S1, S2 on %s with prefix %s."
            % (sourcePath, options.prefix) )
    if options.verbose:
        pass
    if options.destPath is not None:
        pass

if __name__ == '__main__':
    main()
