"""
Reference PyMol scripts
"""

from pymol import cmd
from pymol import stored
from pymol import selector
from chempy import cpv


def colorbyrmsd(mobile, target, doAlign=1, doPretty=1, guide=1, method='super', quiet=0):
    """
http://pymolwiki.org/index.php/ColorByRMSD

Original Authors: Shivender Shandilya; Jason Vertrees
Complete rewrite by Thomas Holder

License: BSD-2-Clause

DESCRIPTION

    Align two structures and show the structural deviations in color to more
    easily see variable regions.

    Colors each mobile/target atom-pair by distance (the name is a bit
    misleading).

    Modifies the B-factor columns in your original structures.

ARGUMENTS

    mobile = string: atom selection for mobile atoms

    target = string: atom selection for target atoms

    doAlign = 0 or 1: Superpose selections before calculating distances
    {default: 1}

    doPretty = 0 or 1: Show nice representation and colors {default: 1}

EXAMPLE

    fetch 1ake 4ake, async=0
    remove chain B
    colorbyrmsd 1ake, 4ake
    """
    doAlign, doPretty = int(doAlign), int(doPretty)
    guide, quiet = int(guide), int(quiet)
    aln, seleboth = '_aln', '_objSelBoth'
    try:
        align = cmd.keyword[method][0]
    except:
        print(' Error: no such method: ' + str(method))
        return
    if guide:
        mobile = '(%s) and guide' % mobile
        target = '(%s) and guide' % target
    try:
        if doAlign:
            # superpose
            align(mobile, target)
        # get alignment without superposing
        align(mobile, target, cycles=0, transform=0, object=aln)
    except:
        print(' Error: Alignment with method %s failed' % (method))
        return
    cmd.select(seleboth, '(%s) or (%s)' % (mobile, target))
    idx2coords = dict()
    cmd.iterate_state(-1, seleboth, 'idx2coords[model,index] = (x,y,z)', space=locals())
    if cmd.count_atoms('?' + aln, 1, 1) == 0:
        # this should ensure that "aln" will be available as selectable object
        cmd.refresh()
    b_dict = dict()
    for col in cmd.get_raw_alignment(aln):
        assert len(col) == 2
        b = cpv.distance(idx2coords[col[0]], idx2coords[col[1]])
        for idx in col:
            b_dict[idx] = b
    cmd.alter(seleboth, 'b = b_dict.get((model, index), -1)', space=locals())
    if doPretty:
        cmd.orient(seleboth)
        cmd.show_as('cartoon', 'byobj ' + seleboth)
        cmd.color('gray', seleboth)
        cmd.spectrum('b', 'blue_red', seleboth + ' and b > -0.5')
    if not quiet:
        print(" ColorByRMSD: Minimum Distance: %.2f" % (min(b_dict.values())))
        print(" ColorByRMSD: Maximum Distance: %.2f" % (max(b_dict.values())))
        print(" ColorByRMSD: Average Distance: %.2f" % (sum(b_dict.values()) / len(b_dict)))
    cmd.delete(aln)
    cmd.delete(seleboth)


cmd.extend('xx_colorbyrmsd', colorbyrmsd)
cmd.auto_arg[0]['xx_colorbyrmsd'] = cmd.auto_arg[0]['align']   # for auto-completion
cmd.auto_arg[1]['xx_colorbyrmsd'] = cmd.auto_arg[1]['align']


def readSymmetry(inFile):
    """
    Link:   https://pymolwiki.org/index.php/BiologicalUnit
    Author: Jason Vertrees <Jason-dot-Vertrees-at-schrodinger_dot_com>, 2010

    Learn:  https://proteopedia.org/wiki/index.php/Biological_unit
        ->  https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies


    This function will read "inFile" and glean the symmetry operations, if any, from it.

    PARAMS
    inFile
        (string) path to PDB file
        
    RETURNS
    matrix
        Array of lists.  One 16-element list per symmetry operation.  Feed this matrix
        into manualSymExp in order to make the other symmetry mates in the biological unit
    """
    # REMARK 350 BIOMTn TAGn X Y Z Tx
    rem, tag, biomt, opcode, x, y, z, tx = list(range(8))
    thePDB = open(inFile, 'rt').readlines()
    matrices = []
    curTrans = -1
    # The transformation is,
    # output = U*input + Tx
    for l in thePDB:
        tokens = l.split()
        if len(tokens)!=8:
            continue
        if tokens[rem]=="REMARK" and tokens[tag]=="350" and tokens[biomt].startswith("BIOMT"):
            if tokens[opcode]!=curTrans:
                # new transformation matrix
                matrices.append([])
                matrices[-1].append( map( lambda s: float(s), tokens[x:]))
            curTrans = tokens[opcode]
    print("Found %s symmetry operators in %s." % (len(matrices), inFile))

    """
    for m in matrices:
        n = cmd.get_unused_name(prefix)
        cmd.create(n, objSel)
        s1 = "%s + (x*%s + y*%s + z*%s)" % (m[0][3], m[0][0], m[0][1], m[0][2])
        s2 = "%s + (x*%s + y*%s + z*%s)" % (m[1][3], m[1][0], m[1][1], m[1][2])
        s3 = "%s + (x*%s + y*%s + z*%s)" % (m[2][3], m[2][0], m[2][1], m[2][2])
        cmd.alter_state(1, n, "(x,y,z) = (%s, %s, %s)" % (s1, s2, s3) )
    """
    return matrices

cmd.extend('xx_readSymmetry', colorbyrmsd)




