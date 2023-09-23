#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import datetime as dtm
import numpy as np
import numpy.ma as ma
import netCDF4 as nc4
import argparse
import logging
import fnmatch
from collections import OrderedDict
import matplotlib as mpl
import matplotlib.pyplot as plt


#--------------------------------------------------
#     l o g g i n g
#
PkgLogger  = logging.getLogger( __name__ )
OUT_HDLR = logging.StreamHandler( sys.stdout )
OUT_HDLR.setFormatter( logging.Formatter( '%(asctime)s %(levelname)s::%(funcName)s:: %(message)s') )
OUT_HDLR.setLevel( logging.INFO )
PkgLogger.addHandler( OUT_HDLR )
PkgLogger.setLevel( logging.INFO )


def pkglog_set_level( level ):
    """Function to set/change logging level of package logging facility

    :level:  log level value, like e.g. logging.INFO, logging.DEBUG, ...
    """
    OUT_HDLR.setLevel(level)
    PkgLogger.setLevel(level)

# ---pkglog_set_level---


def beps_site_results_compare(fname1,fname2,
                              disable_automask=False, stop_on_diff=True, tol=None):
    """Function that compares to BEPS generated NetCDF datafiles
    (that typically contain simulated quantities for one day)

    Returns:
    0: if files have equal data sets
    1: otherwise
    """
    rc = 0
    
    fp1 = nc4.Dataset(fname1)
    fp2 = nc4.Dataset(fname2)

    #-- disable automatic scaling/masking (at least until missing_value issue is solved)
    if disable_automask:
        fp1.set_auto_mask(False)
        fp2.set_auto_mask(False)

    #-- ensure both NetCDF files contain the same datasets
    vlist1 = list(fp1.variables.keys())
    vlist2 = list(fp2.variables.keys())
    if not vlist1==vlist2:
        msg = "incompatible variable lists, cannot continue!"
        msg += " (vlist1 +++{}+++ vlist2 +++{}+++".format(vlist1, vlist2)
        PkgLogger.fatal(msg)
        raise RuntimeError(msg)
    else:
        #-- convert time-point values to dates
        r = fp1.variables['time'].units.split()
        ref_time = None
        for fmt in ['%Y%m%dT%H:%M:%S', '%Y-%m-%dT%H:%M:%S']:
            try:
                ref_time = dtm.datetime.strptime(r[-1], fmt)
                break
            except ValueError:
                pass
        if ref_time==None:
            msg = "reference time could not be parsed correctly from string ***{}***".format(r)
            PkgLogger.fatal(msg)
            raise RuntimeError(msg)

        r[-1] = ref_time.strftime('%Y-%m-%dT%H:%M:%S')
        tunit = ' '.join(r)
        dates_a = nc4.num2date(fp1.variables['time'][:], tunit)

    #-- simulated day
    act_day = dates_a[0].strftime('%Y%m%d')
    #-- loop over variables
    reldiff_dct = OrderedDict()
    for k in vlist1:
        if k in ['time','nsite']:
            assert np.all(fp1.variables[k][:]==fp2.variables[k][:]), \
                "inconsistency at aux variable +++{}+++".format(k)
            continue
        msg = "start comparison for dataset +++{}+++...".format(k)
        PkgLogger.info(msg)
        #-- datasets as numpy arrays
        var1 = fp1.variables[k][:]
        var2 = fp2.variables[k][:]
        #-- we do *NOT* expect fill-values in BEPS generated files
        if not disable_automask:
            if ma.count_masked(var1)>0:
                msg = "file={} var={} with unexpected fill-values ({})".format(
                    fname1, k, ma.count_masked(var1))
                PkgLogger.error(msg)
                rc = 1
                if stop_on_diff:
                    break
            elif ma.count_masked(var2)>0:
                msg = "file={} var={} with unexpected fill-values ({})".format(
                    fname2, k, ma.count_masked(var2))
                PkgLogger.error(msg)
                rc = 1
                if stop_on_diff:
                    break
        #-- absolute difference
        rc = 0
        var_diff = ma.abs(var1-var2)
        if var_diff.max()==0.:
            reldiff = 0.
        else:
            v1 = np.abs(var1).ravel()
            v2 = np.abs(var2).ravel()
            vdif = var_diff.ravel()
            vmx = np.maximum(v1,v2)
            cnd = (vmx>0)
            if len(cnd)>0:
                reldiff = np.max(vdif[cnd]/vmx[cnd])
            else:
                reldiff = 0.
        reldiff_dct[k] = reldiff
        #-- check absolute difference
        if var_diff.max()>0.:
            imx = var_diff.argmax()
            msg = "detected non-zero differences for +++{}+++ (min/mean/max = {}/{}/{})".format(
                k, var_diff.min(), var_diff.max(), var_diff.mean())
            PkgLogger.error(msg)
            rc = 1
            if stop_on_diff:
                break
        else:
            msg = "dataset={}, differences min/max/mean={}/{}/{}".format(
                k, var_diff.min(), var_diff.max(), var_diff.mean())
            PkgLogger.debug(msg)
        #-- logging
        msg = "...comparison DONE (+++{}+++, rc={})".format(k,rc)
        PkgLogger.info(msg)

    #-- close handles
    fp1.close()
    fp2.close()

    return (rc,reldiff_dct)
# ---end-of-beps_site_results_compare


def subcmd_cmp_result_dirs(options):
    """Function ...
    """
    #-- grab options
    dir1,dir2 = options.dirs
    disable_automask = options.disable_automask
    stop_on_diff = options.stop_on_diff

    #-- logging
    msg = "comparison direcories ***{}*** ***{}***".format(dir1,dir2)
    PkgLogger.info(msg)

    #-- BEPS file name pattern (dayly output files expected)
    ptn = 'beps_site_????????.nc'
    dir1_flst = sorted(fnmatch.filter(os.listdir(dir1),ptn))
    dir2_flst = sorted(fnmatch.filter(os.listdir(dir2),ptn))

    #-- check both directories contain the same NetCDF files
    if dir1_flst!=dir2_flst:
        msg = "inconsistent list of NetCDF result files (dir1 ***{}***,  dir2 ***{}***)".format(
            dir1_flst, dir2_flst)
        PkgLogger.fatal(msg)
        sys.exit(1)

    #-- comparison-loop
    reldiff_dct = None
    for i,f in enumerate(dir1_flst):
        fqf1 = os.path.join(dir1,f)
        fqf2 = os.path.join(dir2,f)
        msg = "START comparison of files ***{}*** ***{}***...".format(fqf1, fqf2)
        PkgLogger.info(msg)
        rc,act_rdiff_dct = beps_site_results_compare(fqf1, fqf2,
                                                     disable_automask=disable_automask,
                                                     stop_on_diff=stop_on_diff)
        if i==0:
            reldiff_dct = act_rdiff_dct
        else:
            for k in reldiff_dct.keys():
                rd1 = reldiff_dct[k]
                rd2 = act_rdiff_dct[k]
                reldiff_dct[k] = max(rd1, rd2)
        #-- stop on first difference detected
        if rc!=0:
            msg = "detected difference at files ***{}***   ***{}***.".format(
                fqf1, fqf2)
            PkgLogger.error(msg)
            if stop_on_diff:
                msg = "Stop further comparison!"
                PkgLogger.info(msg)
                sys.exit(rc)
                break
        msg = "...comparison of all files DONE"
        PkgLogger.info(msg)

    #-- maximal relative difference
    if options.outname!=None:
        outname = options.outname
    else:
        outname = 'reldiff_max.txt'
    with open(outname,'w') as fp:
        fp.write("{:<20s} {:<22s}".format("dataset", "max relative difference") + '\n')
        for k,v in reldiff_dct.items():
            fp.write("{:<20s} {:22.15e}".format(k,v) + '\n')
        msg = "written file ***{}***".format(outname)
        PkgLogger.info(msg)
# ---end-of-subcmd_cmp_result_dirs


def create_argument_parser(progname=None):
    """Function to create the command line parser

    :progname: file/programme name of script

    :return    argument parser suitable for use in main routine
    :rtype:    ArgumentParser instance
    """
    from argparse import ArgumentParser

    def _add_io_options(aparser):
        aparser.add_argument( '--outname',
                              help="""write output to this file""" )
        # aparser.add_argument( '--outdir',
        #                       help="""where to place generated file(s).""" )
    # ---_add_io_options---

    def _add_plot_options(aparser):
        aparser.add_argument( '--dpi',
                              type=int,
                              default=150,
                              help="""dots-per-inch in generated plot(s)""" )
        aparser.add_argument( '--title',
                              help="""title at top of generated plot (Note: maybe augmented in case multiple plots need to be generated.)""" )
    # ---_add_plot_options---

    #-- create main instance
    parser = ArgumentParser( prog=progname, usage=globals()['__doc__'] )

    #-------------------
    #     c o m m o n   o p t i o n s
    #
    parser.add_argument( '--verbose','-v',
                         action='store_true',
                         help="""dump some more debugging messages to terminal.""" )

    #----------------------------
    #     s u b c o m m a n d s
    #
    subparsers = parser.add_subparsers( title='Available Subcommands',
                                        metavar='CMDS',
                                        description='',
                                        dest='subcmds',
                                        help=''
                                        )


    #-------------------
    #      c m p _ r e s u l t _ d i r s
    #
    xparser = subparsers.add_parser( 'cmp_result_dirs',
                                     help="""compare two directories resulting from BEPS site-level run. Assuming that all existing NetCDF datafiles contain the same datasets, the maximal relative diffference per dataset (over all time-points) is written to file. """ )
    xparser.add_argument('dirs',
                         nargs=2,
                         help="""(fully qualified) names of both resulting directories""")
    xparser.add_argument( '--disable_automask',
                          action='store_true',
                          help="""whether to disable automatic masking of fill-values when reading from NetCDF file (enable this option to address 'missing_value' issue in BEPS generated files.""" )
    xparser.add_argument('--no-stop-on-diff',
                         dest='stop_on_diff',
                         action='store_false',
                         help="""whether to continue comparison of directory/files when differences have been detected.""")
    _add_io_options(xparser)


    
    return parser
# ---end-of-create_argument_parser


def main(options):

    tstart = dtm.datetime.now()

    #-------------------
    #    c m p _ r e s u l t _ d i r s
    #
    if options.subcmds=='cmp_result_dirs':
        subcmd_cmp_result_dirs(options)

    #-------------------
    #    c m p _ j a c o b i a n s
    #
    if options.subcmds=='cmp_jacobians':
        subcmd_cmp_jacobians(options)


    #-- stop timing
    tend = dtm.datetime.now()
    elapsed = (tend-tstart).total_seconds()

    #-- logging
    msg = "subcommand '{}' terminated (elapsed={:.1f}[s])".format(options.subcmds, elapsed)
    PkgLogger.info(msg)
# ---end-of-main


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#                    M A I N
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':
    import datetime as dtm

    progname = os.path.basename(__file__) #determine filename

    #-----------------------------
    #          P R O G R A M   S T A R T
    #
    fmt = "%Y-%m-%dT%H:%M:%S.%f"
    ttstart = dtm.datetime.now()
    PkgLogger.info("{}::PROGRAM START::{}".format(progname, ttstart.strftime(fmt)))
    PkgLogger.info( "command-line: {}".format(' '.join(sys.argv)) )

    #          p a r s e   c o m m a n d   l i n e
    #
    parser = create_argument_parser(progname=progname)
    options = parser.parse_args()
    if options.verbose:
        pkglog_set_level(logging.DEBUG)

    #-- call main method
    main(options)

# ---M A I N---
