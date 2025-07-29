"""MDAnalysis

Unit:   https://userguide.mdanalysis.org/stable/units.html

    lenght      Angstrom
    time        ps
    energy      kJ/mol
    charge      elementary
    mass        atomic
    angle       degree


Info:   https://userguide.mdanalysis.org/stable/topology_system.html#topology-attributes

    Be aware that `id` can be repeated if number of atoms are bigger than 99999, use `ix` instead,
    other cannoical attribute are (index,indices), (resindex,resindices), (segindex,segindices)

    AtomGroup = u.select_atoms([CMD | Slice])

    # iteration is on timestamp, coordinates follow
    u.trajectory[10].frame      # 10

    # selections are defined on static properties that do not change
    # however, some selections contain distance-dependent queries, `updating=True` needed
    dynamic = u.select_atoms('around 2 resname ALA', updating=True)
        # <class 'MDAnalysis.core.groups.UpdatingAtomGroup'>
    static = u.select_atoms('around 2 resname ALA')
        # <class 'MDAnalysis.core.groups.AtomGroup'>

    # write to file
    ca = u.select_atoms('name CA')
    ca.write('calphas.gro')

    ca = u.select_atoms('name CA')
    with mda.Writer('calphas.xtc', ca.n_atoms) as w:
        for ts in u.trajectory:
            w.write(ca)

    # if not provided, atomic mass is guessed, thus check them before analysis
    u.atoms.groupby('masses')
    u.atoms.groupby('types')



"""


import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
from numpy.linalg import norm, svd
import tqdm

import os
import json
import shutil
import hashlib


PROTEIN_RING_ATOMS = {
    'TRP': [['CD2', 'CE3', 'CZ3', 'CH2', 'CZ2', 'CE2'],
            ['CD2', 'CE2', 'NE1', 'CD1', 'CG']],        # special for 5-atom ring
    'PHE': ['CD1', 'CE1', 'CZ', 'CE2', 'CD2', 'CG'],
    'TYP': ['CE1', 'CZ', 'CE2', 'CD2', 'CG', 'CD1'],
    'HIS': ['ND1', 'CE1', 'NE2', 'CD2', 'CG'],
}
PROTEIN_RING_ATOMS['HIP'] = PROTEIN_RING_ATOMS['HID'] = PROTEIN_RING_ATOMS['HIS']


def allclose(a, b, tol=None) -> bool:
    """recursively compare a and b on given tolerance"""
    tol = tol if tol else 0.000001
    if a is None and b is None:
        return True
    if isinstance(a, bool) and isinstance(b, bool):
        return a == b
    if isinstance(a, (int,float)) and isinstance(b, (int,float)):
        if abs(a-b) > tol:
            return False
        return True
    if isinstance(a, (tuple,list)) and isinstance(b, (tuple,list)):
        if len(a) != len(b):
            return False
        if not a:           # for empty object
            return True
        for vi,vj in zip(a,b):
            if not allclose(vi, vj, tol):
                return False
        return True
    if isinstance(a, dict) and isinstance(b, dict):
        if set(list(a.keys())) != set(list(b.keys())):
            return False
        for k in a.keys():
            if not allclose(a[k], b[k], tol):
                return False
        return True
    return False


def xhash(obj) -> int:
    def serialize(o):
        if isinstance(o, (int, float, bool, type(None))):
            return repr(o)
        elif isinstance(o, str):
            return f'str:{o}'
        elif isinstance(o, (list, tuple)):
            return f"{type(o).__name__}:[{','.join(serialize(i) for i in o)}]"
        elif isinstance(o, dict):
            # sort keys to ensure determinism
            items = sorted((serialize(k), serialize(v)) for k, v in o.items())
            return f"dict:{{{','.join(f'{k}:{v}' for k, v in items)}}}"
        else:
            raise TypeError(f"Unsupported type: {type(o)}")

    serialized = serialize(obj)
    serialbytes = hashlib.sha256(serialized.encode('utf-8')).digest()
    return int.from_bytes(serialbytes[:8], byteorder='big', signed=False)


class PiPiStack:
    """
    Args:
        u (mda.Universe):
        groupA | groupB (dict): {resname: [atname1, atname2, ...]}  : ring atom names
        groupA_add_protein_residues | groupB_add_protein_residues (bool):
            whether added standard protein residues to each group

        groups (dict): self rings pi-pi stacking, highest priority

        bak_stepsize (int): calculate every this stepsize
        bak_stepsave (int): bak result every this step, `filename-pistack.json`
        bak_continue (bool): whether continue run to save time

        cutoff_distance (float): distance between two rings' center
        cutoff_angle (float): acute angle between two rings' normal vectors
        cutoff_offset (float): distance for centroid points projection on each rings' plane

    Examples:
        # for non-standard ring atoms, we have/define,
        >>> non_std_ring_atoms = {
                'RS1':  ['a1', 'a2', 'a3', 'a4', 'a5'],
                'RS2':  [['b1', 'b2', 'b3', 'b4', 'b5', 'b6'],
                         ['c1', 'c2', 'c3', 'c4', 'c5']]
            }

        # case 1, calculate all pi-pi stacking for all rings
        >>> ps = PiPiStack(
                u, groups=non_std_ring_atoms,
                groupA_add_protein_residues=True, groupB_add_protein_residues=True
            )

        # case 2, calculate only non-standard proteins pi-pi stacking
        >>> ps = PiPiStack(
                u, groups=non_std_ring_atoms,
                groupA_add_protein_residues=False, groupB_add_protein_residues=False
            )

        # case 3, calculate pi-pi stacking between non-standard protein and all rings
        >>> ps = PiPiStack(
                u, groups=non_std_ring_atoms,
                groupA_add_protein_residues=False, groupB_add_protein_residues=True
            )

        # case 4, calculate pi-pi stacking between non-standard protein and standard proteins
        >>> ps = PiPiStack(
                u, groupA=non_std_ring_atoms,
                groupA_add_protein_residues=False, groupB_add_protein_residues=True
            )

    Checkpoint:
        JSON file, keys:
            bak_file | bak_prevfile : checkfile name
            bak_stepsave : frequency to be backed up
            bak_start : the last step to be processing
            bak_ix (List[List[int]]) : ring atoms ix
            '0', '1', ..., 'n': correspond frame index
                [ [ringA_ix, ringB_ix, distance, pi_pi_type, angle, offsetA, offsetB, plane_angle] ]
                where, pi_pi_type: 0: parallel, 1: T-shaped, 2: transition
    """
    def __init__(
        self, u, groupA=None, groupB=None, groups=None, bak_stepsize=None,
        bak_stepsave=None, bak_continue=True, bak_check_ring_planarity=None,
        groupA_add_protein_residues=True, groupB_add_protein_residues=True,
        cutoff_distance=5.5, cutoff_angle=30, cutoff_offset=3.5
    ):
        self.u = u
        self.groupA = groupA if groupA else {}
        self.groupB = groupB if groupB else {}
        if groups:
            if self.groupA:
                print('PiPiStack: overwrite groupA by groups')
            if self.groupB:
                print('PiPiStack: overwrite groupB by groups')
            self.groupA = groups
            self.groupB = groups
        self.cutoff_distance = cutoff_distance
        self.cutoff_angle = cutoff_angle
        self.cutoff_offset = cutoff_offset
        self.groupA_add_protein_residues = True if groupA_add_protein_residues is True else False
        self.groupB_add_protein_residues = True if groupB_add_protein_residues is True else False
        self._ring_groups = None
        self.bak_stepsave = bak_stepsave
        self.bak_continue = True if bak_continue is True else False
        self.bak_start = 0
        self.bak_stepsize = bak_stepsize if bak_stepsize else 1
        self.bak_check_ring_planarity = True if bak_check_ring_planarity is True else False
        self.final = {}

    def get_ring_atoms_ix(self, u, ring_atoms=None) -> list:
        if not ring_atoms: ring_atoms = PROTEIN_RING_ATOMS
        ring_atoms_ix = []
        for r in u.residues:
            if r.resname in ring_atoms:
                p = ring_atoms[r.resname]
                if isinstance(p[0],str):
                    p = [p, ]
                for s in p:
                    ats = ' '.join(s)
                    g = u.select_atoms('resindex {:} and name {:}'.format(r.resindex, ats))
                    ring_atoms_ix.append(g.ix.tolist())     # list[int] type!!!
        return ring_atoms_ix

    def _get_pi_pi_ix_groups(self):
        if self.groupA_add_protein_residues:
            ring_atoms_A = {**self.groupA, **PROTEIN_RING_ATOMS}
        else:
            ring_atoms_A = self.groupA
        ring_ix_A = self.get_ring_atoms_ix(self.u, ring_atoms_A)

        if self.groupB_add_protein_residues:
            ring_atoms_B = {**self.groupB, **PROTEIN_RING_ATOMS}
        else:
            ring_atoms_B = self.groupB
        ring_ix_B = self.get_ring_atoms_ix(self.u, ring_atoms_B)

        print('PiPiStack: get ring ix')
        groups = []
        for ia in tqdm.tqdm(ring_ix_A):
            sa = set(ia)
            for ib in ring_ix_B:
                if len(sa) != len(sa.union(ib)):    # remove self repeats
                    groups.append([sorted(ia),sorted(ib)])

        """
        # too time consuming
        keep = []
        done = []
        print('PiPiStack: remove ring ix duplicates')
        for g in tqdm.tqdm(groups):
            p1 = g[0] + g[1]
            p2 = g[1] + g[0]
            if p1 not in done and p2 not in done:   # remove duplicates, (ia,ib) === (ib,ia)
                keep.append(g)
                done.append(p1)
                done.append(p2)
        """

        keep = []
        done = set()
        print('PiPiStack: remove ring ix duplicates')
        for g in tqdm.tqdm(groups):
            p1 = xhash(g[0]+g[1])
            p2 = xhash(g[1]+g[0])
            if p1 not in done and p2 not in done:   # remove duplicates, (ia,ib) === (ib,ia)
                keep.append(g)
                done.add(p1)
                done.add(p2)
        print('  => number of pi ring groups: {:}'.format(len(keep)))

        return keep

    def _single_frame(self):
        # format:
        #   List[ [ringA_ix, ringB_ix, distance, pi_pi_type, angle, offsetA, offsetB, plane_angle] ]
        #   where, pi_pi_type: 0: parallel, 1: T-shaped, 2: transition
        keep = []
        if not self._ring_groups:
            self._ring_groups = self._get_pi_pi_ix_groups()

        boxfull = self.u.dimensions[:3]
        boxhalf = boxfull / 2
        for (ia,ib) in self._ring_groups:
            ra = self.u.atoms[ia].positions     # now, process PBC
            ca = ra.mean(axis=0)
            da = ra - ca
            x1, y1 = np.where(da > boxhalf)
            ee = False
            if len(x1):
                ra[x1, y1] -= boxfull[y1]
                ee = True
            x2, y2 = np.where(da < -boxhalf)    # minus
            if len(x2):
                ra[x2, y2] += boxfull[y2]
                ee = True
            if ee:
                ca = ra.mean(axis=0)

            rb = self.u.atoms[ib].positions
            cb = rb.mean(axis=0)
            db = rb - cb
            x1, y1 = np.where(db > boxhalf)
            ee = False
            if len(x1):
                rb[x1, y1] -= boxfull[y1]
                ee = True
            x2, y2 = np.where(db < -boxhalf)    # minus
            if len(x2):
                rb[x2, y2] += boxfull[y2]
                ee = True
            if ee:
                cb = rb.mean(axis=0)

            vab = cb - ca
            pd = norm(vab)
            if pd > self.cutoff_distance: continue

            auv, asv, avh = svd(ra - ca)
            na = avh[2]
            buv, bsv, bvh = svd(rb - cb)
            nb = bvh[2]

            angle = np.arccos( np.dot(na, nb) ) / np.pi * 180.0
            angle = min(angle, 180-angle)
            if angle <= self.cutoff_angle:
                kp = [ia,ib,pd,0,angle,0,0,0]
            elif angle >= 90 - self.cutoff_angle:
                kp = [ia,ib,pd,1,angle,0,0,0]
            else:
                continue

            ta = np.dot(vab, na)
            oa = norm( vab - ta * na )
            if oa > self.cutoff_offset: continue

            tb = np.dot(-vab, nb)
            ob = norm( -vab - tb * nb )
            if ob > self.cutoff_offset: continue
            kp[5] = oa
            kp[6] = ob

            av = np.arccos( ta / pd ) / np.pi * 180.0
            bv = np.arccos( tb / pd ) / np.pi * 180.0
            av = min(av, 180-av)
            bv = min(bv, 180-bv)
            if kp[3] == 0 and av <= 30 and bv <= 30:
                pass
            elif kp[3] == 1 and ((av >= 60 and bv <= 30) or (av <= 30 and bv >= 60)):
                pass
            else:
                kp[3] = 2
            kp[-1] = (av+bv)/2
            keep.append(kp)

        return keep

    def check_ring_planarity(self, xyzs, tolerance=0.3):
        coords = np.array(xyzs)
        centroid = coords.mean(axis=0)
        uv, sv, vh = svd(coords - centroid)
        if sv[2] > tolerance:
            return False
        return True

    def run(self):
        final = self._prepare()
        bak_file = final['bak_file']
        bak_prevfile = final['bak_prevfile']
        final['bak_stepsave'] = self.bak_stepsave
        final['bak_start'] = self.bak_start
        print('PiPiStack: stepsize {:}, bak_start {:}, bak_stepsize {:}'.format(
            self.bak_stepsave, self.bak_start, self.bak_stepsize)
        )

        if not self._ring_groups:
            self._ring_groups = self._get_pi_pi_ix_groups()
        if not self._ring_groups:
            print('PiPiStack: no ring found')
            return

        if self.bak_check_ring_planarity:
            print('PiPiStack: check rings planarity')
            done = set()
            for (ia,ib) in self._ring_groups:
                ha = xhash(ia)
                if ha not in done:
                    done.add(ha)
                    ra = self.u.atoms[ia].positions
                    if not self.check_ring_planarity(ra):
                        print('  -> not-in-plane: {:}'.format(ia))
                hb = xhash(ib)
                if hb not in done:
                    done.add(hb)
                    rb = self.u.atoms[ib].positions
                    if not self.check_ring_planarity(rb):
                        print('  -> not-in-plane: {:}'.format(ib))
            print()

        print('PiPiStack: calculating...')
        if 'bak_ix' in final:
            bak_ix = {xhash(v):i for i,v in enumerate(final['bak_ix'])}
        else:
            final['bak_ix'] = []
            bak_ix = {}
        c = 0
        for ts in tqdm.tqdm(self.u.trajectory[self.bak_start::self.bak_stepsize]):
            keep = self._single_frame()
            if not keep: continue

            c += 1
            for g in keep:      # trim, convert np.float_ to float!!!
                g[2] = round(float(g[2]), 3)
                g[4] = round(float(g[4]), 3)
                g[5] = round(float(g[5]), 3)
                g[6] = round(float(g[6]), 3)
                g[7] = round(float(g[7]), 3)

                ha = xhash(g[0])
                if ha in bak_ix:
                    g[0] = bak_ix[ha]                               # use exist
                else:
                    final['bak_ix'].append(g[0])
                    bak_ix[ha] = g[0] = len(final['bak_ix']) - 1    # update

                hb = xhash(g[1])
                if hb in bak_ix:
                    g[1] = bak_ix[hb]
                else:
                    final['bak_ix'].append(g[1])
                    bak_ix[hb] = g[1] = len(final['bak_ix']) - 1

            x = int(self.u.trajectory.frame)
            final[x] = keep
            if c % self.bak_stepsave == 0:
                if os.path.isfile(bak_file):
                    shutil.move(bak_file, bak_prevfile)

                final['bak_start'] = x
                with open(bak_file, 'w') as f:
                    json.dump(final, f)

        if len(self.u.trajectory) > final['bak_start'] + 1:      # left
            final['bak_start'] = len(self.u.trajectory)
            with open(bak_file, 'w') as f:
                json.dump(final, f)
        self.final = final

    def _prepare(self):
        base, ext = os.path.splitext(self.u.filename)
        bak_file = base + '-pistack.json'
        bak_prevfile = base + '-pistack_prev.json'
        finalbak = {}
        if os.path.isfile(bak_file):
            with open(bak_file, 'r') as f:
                finalbak = json.load(f)

        if not self.bak_stepsave:
            if 'bak_stepsave' in finalbak:
                self.bak_stepsave = finalbak['bak_stepsave']
            else:
                self.bak_stepsave = int(5000000 / len(self.u.atoms))
                self.bak_stepsave = min(15, self.bak_stepsave)

        if 'bak_start' in finalbak:
            self.bak_start = finalbak['bak_start']
        
        if 'bak_stepsize' in finalbak:
            self.bak_stepsize = finalbak['bak_stepsize']
        else:
            finalbak['bak_stepsize'] = self.bak_stepsize

        finalbak.update({'bak_file': bak_file, 'bak_prevfile': bak_prevfile})
        return finalbak

    def write_pi_pi_stack_to_xyzfile(self, frame=None, ptype=None):
        """
        args:
            frame (None | int | list[int]): which frame to be wrote, starts from 0
            ptype: None | 'all' | 'p' (parallel) | 't' (T-shaped) | 's' (transition)

        explain:
            (frame=None, ptype=None)  :  write all frames all types pi-pi stacking
            (frame=None, ptype='p')   :  write all frames parallel pi-pi stacking
            (frame=[1,4],ptype='t')   :  write 2nd and 5th frames T-shaped pi-pi stacking
            (frame=None, ptype=['p','t'])   :  write all frames parallel and T-shaped pi-pi stacking
        """
        if not self.final:
            self.final = self._prepare()
        if not self.final:
            print('PiPiStack: no results')
            return

        if frame:
            if isinstance(frame, (int,str)):
                frame = [str(frame), ]
        else:
            frame = [str(i) for i in range(len(self.u.trajectory))]

        keys = set(frame).intersection(list(self.final.keys()))
        if not keys:
            print('PiPiStack: no selections')
        keys = sorted(list(map(int,keys)))

        if not ptype:
            sel = set([0, 1, 2])
        else:
            if isinstance(ptype, str): ptype = [ptype, ]
            sel = set()
            for s in ptype:
                if s.lower() in ['a', 'all']:
                    sel = [0, 1, 2]
                    break
                elif s.lower() in ['p', 'parallel']:
                    sel.add(0)
                elif s.lower() in ['t', 't-shaped', 't-shape', 'tshaped', 'tshape']:
                    sel.add(1)
                else:
                    sel.add(2)
            if not sel: sel = set([0, 1, 2])
        lp = ['Parallel', 'T-shaped', 'Transition', 'All']
        if len(sel) == 3:
            selp = lp[-1]
        else:
            selp = '|'.join([lp[i] for i in sorted(sel)])

        base, ext = os.path.splitext(self.final['bak_file'])
        xyzfile = base + '.xyz'
        print('PiPiStack: writing to file: {:}'.format(xyzfile))
        bak_ix = self.final['bak_ix']

        with open(xyzfile, 'wt') as f:
            for k in keys:
                frame = self.u.trajectory[k]
                boxfull = self.u.dimensions[:3]
                boxhalf = boxfull / 2
                for i,g in enumerate(self.final[str(k)]):
                    if g[3] not in sel: continue

                    ia = bak_ix[g[0]]
                    ib = bak_ix[g[1]]
                    xx = [*ia, *ib]
                    f.write('{:}\n'.format(len(xx)))    # number of atoms
                    f.write('f{:}-i{:}-{:},d-{:},a-{:},oa-{:},ob-{:},plane-{:}\n'.format(
                        str(k),i,selp,g[2],g[4],g[5],g[6],g[7]
                    ))

                    gt = self.u.atoms[xx].types
                    gp = self.u.atoms[xx].positions     # now, process PBC
                    dd = gp - gp.mean(axis=0)
                    x1, y1 = np.where(dd > boxhalf)
                    if len(x1):
                        gp[x1, y1] -= boxfull[y1]
                    x2, y2 = np.where(dd < -boxhalf)    # minus
                    if len(x2):
                        gp[x2, y2] += boxfull[y2]

                    for t,p in zip(gt, gp):
                        f.write('{:3}     {:8.3f}  {:8.3f}  {:8.3f}\n'.format(t, *p))


if __name__ == '__main__':

    rings = {
        'PPE': [
            ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
            ['C9', 'C10', 'C11', 'C12', 'C13', 'C14'],
            ['C15', 'C17', 'C18', 'C19', 'C20', 'C21'],
            ['C16', 'C22', 'C23', 'C24', 'C25', 'C26'],
            ['C29', 'C30', 'C31', 'N32', 'C33', 'C34'],
            ['N39', 'N40', 'N41', 'C42', 'C43']
        ]
    }

    print('starting->\n')
    U = mda.Universe('pp4-min.gro', 'pp4-prod-full-100ns.xtc')
    print('processing->\n')
    ps = PiPiStack(U, groups=rings, bak_check_ring_planarity=True, bak_stepsize=100)
    #ps.run()


    ps.write_pi_pi_stack_to_xyzfile(frame=100,ptype=None)



    print('DONE')




