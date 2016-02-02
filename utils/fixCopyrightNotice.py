#! /usr/bin/env python

"""(Ostensibly) fix copyright notices in files.

Actually, this sript will simply replace a block of text in a file from one
string to another.  It will only do this once though, i.e. not globally
throughout the file.  It writes a backup file and then does an os.rename()
dance for atomicity.

Usage: fixnotices.py [options] [filenames]
Options:
    -h / --help
        Print this message and exit

    --newnotice=file
        Use the notice in the file as the new (replacement) string, instead of
        the hard coded value in the script.

    --dry-run
        Don't actually make the changes, but print out the list of files that
        would change.  When used with -v, a status will be printed for every
        file.

    -v / --verbose
        Print a message for every file looked at, indicating whether the file
        is changed or not.
"""

import os
import sys
import getopt
import re

NEW_NOTICE = """/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
"""

DRYRUN = 0
VERBOSE = 0

#regexp = re.compile('/\*(.|[\r\n])*\*/')
regexp = re.compile('/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/')

def usage(code, msg=''):
    print __doc__ % globals()
    if msg:
        print msg
    sys.exit(code)


def main():
    global DRYRUN, NEW_NOTICE, VERBOSE
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hv',
                                   ['help', 'oldnotice=', 'newnotice=',
                                    'dry-run', 'verbose'])
    except getopt.error, msg:
        usage(1, msg)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage(0)
        elif opt in ('-v', '--verbose'):
            VERBOSE = 1
        elif opt == '--dry-run':
            DRYRUN = 1
        elif opt == '--oldnotice':
            fp = open(arg)
            OLD_NOTICE = fp.read()
            fp.close()
        elif opt == '--newnotice':
            fp = open(arg)
            NEW_NOTICE = fp.read()
            fp.close()

    for arg in args:
        process(arg)

def process(file):
    f = open(file)
    data = f.read()
    f.close()
    found = regexp.match(data)
    if found is None:
        if VERBOSE:
            print 'no change:', file
        return
    elif DRYRUN:
        found = found.group(0)
        print '   change:', file
        if VERBOSE:
            print '      old:', found
            print '      new:', NEW_NOTICE
    if DRYRUN:
        # Don't actually change the file
        return
    data = regexp.sub(NEW_NOTICE,data,count=1)
    new = file + ".new"
    backup = file + ".bak"
    f = open(new, "w")
    f.write(data)
    f.close()
    os.rename(file, backup)
    os.rename(new, file)


if __name__ == '__main__':
    main()
