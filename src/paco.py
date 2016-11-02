#!/usr/bin/env python
import optparse

def main():
    p = optparse.OptionParser()
    p.add_option('--sourcePath', '-s', default="..\\run1")
    p.add_option('--destPath', '-d', default="..\\run1\\data")
    options, arguments = p.parse_args()
    print 'Mangle from %s to %s' % (options.sourcePath, options.destPath)


if __name__ == '__main__':
    main()
