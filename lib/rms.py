from numpy import *

def rms_fit(ref_c,c,w=None):
    if ( w != None ):
        sum_w = sum(w)

        # move geometric center to the origin
        ref_trans = sum( ref_c * w[:,newaxis], axis=0 ) / sum_w
        ref_c = ref_c - ref_trans
        c_trans = sum( c*w[:,newaxis], axis=0 ) / sum_w
        c = c - c_trans

        # covariance matrix
        C = dot( (c*w[:,newaxis]).T, ref_c ) / ( sum_w/w.size )
    else:
        # move geometric center to the origin
        ref_trans = average( ref_c, axis=0 )
        ref_c = ref_c - ref_trans
        c_trans = average( c, axis=0 )
        c = c - c_trans

        # covariance matrix
        C = dot( c.T, ref_c )

    # Singular Value Decomposition
    ( r1, s, r2 ) = linalg.svd(C)

    # compute sign (remove mirroring)
    if ( linalg.det(C) < 0 ):
        r2[2,:] *= -1.0
    U = dot( r1, r2 )

    return ( c_trans, U, ref_trans )

def rmsd(c1,c2,w=None):
    ( c_trans, U, ref_trans ) = rms_fit( c1, c2, w )
    new_c2 = dot( c2-c_trans, U ) + ref_trans

    rmsd = 0.0
    if ( w != None ):
        rmsd = sqrt( sum( sum( ( c1 - new_c2 )**2, axis=1 ) * w ) / sum(w) )
    else:
        rmsd = sqrt( average( sum( ( c1 - new_c2 )**2, axis=1 ) ) )

    return rmsd, new_c2
