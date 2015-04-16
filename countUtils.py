#-------------------------------------------------------------------------------

@profile
def calculateDnByDs(bins,counts,eucl=False,verbose=False,idl_style=True,
                    errors=False,bright=False,bright_errors=False,
                    return_all=False):
    """

    This function expects bins to be in uJy
    The output units are uJy^-1, or uJy^1.5 if eucl=True

    The function is area-agnostic: sr^-1 -> sr^-1, deg^-2 -> deg^-2 etc.

    Be careful with units:
       e.g. if bins are in uJy, dn_by_ds will be uJy^-1   = 10^6 Jy
                                dn_ny_ds_eucl in uJy^1.5  = 10^-9 Jy
    etc.
    """

    dn=counts
    Smed=medianArray(bins)
    ds=numpy.absolute(numpy.gradient(Smed))

    if idl_style:
        # i.e. the way DERIV.PRO does it
        ds[0] = abs(-3.0*Smed[0] + 4.0*Smed[1] - Smed[2] / 2.0)
        ds[-1] = abs(3.0*Smed[-1] - 4.0*Smed[-2] + Smed[-3] / 2.0)

    if verbose:
        print 'dn',dn
        print 'ds',ds

    if return_all:
        return dn/ds, (Smed**2.5) * dn/ds,\
          (Smed**2.5) * numpy.sqrt(dn)/ds,\
          (Smed**2.0) * dn/ds,\
          (Smed**2.0) * numpy.sqrt(dn)/ds

    if eucl:
        dn_by_ds = (Smed**2.5) * dn/ds
    elif errors:
        dn_by_ds = (Smed**2.5) * numpy.sqrt(dn)/ds
    elif bright:
        dn_by_ds = (Smed**2.0) * dn/ds
    elif bright_errors:
        dn_by_ds = (Smed**2.0) * numpy.sqrt(dn)/ds
    else:
        dn_by_ds = dn/ds

    return dn_by_ds

#-------------------------------------------------------------------------------
