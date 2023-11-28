#!/usr/bin/env python3

from rdkit import Chem
#from rdkit.Chem import Draw

from abc import ABC, abstractmethod
from collections import OrderedDict
import itertools
import pathlib
import subprocess
from typing import Iterable
import weakref
import os
import io


__doc__ = """XZ GMXTopParser based on https://github.com/janjoswig/MDParser/tree/main"""

__all__ = [
    'SUPPORTED_DIRECTIVES',
    'GromacsTopParser',
    'GromacsTopOps',
    'FEATURES',
    '__doc__',
    '__version__',
]

FEATURES = [
    'version 0.1.0  : XZ GMXTopParser',
    'version 0.2.0  : Add more classes',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION


def _trim_locals(d):
    return {k:v for k, v in d.items() if k not in ("self", "__class__")}


def _make_formatter(f):
    def formatter(value):
        return f"{value:{f}}"
    return formatter


class NodeValue(ABC):
    """Abstract base class for node value types"""
    _count = 0
    _node_key_name = "abstract"

    @classmethod
    def reset_count(cls, value: int = None):
        if value is None: value = 0
        cls._count = value

    @classmethod
    def increase_count(cls):
        cls._count += 1

    @property
    def count(self):
        return self._count

    def __init__(self):
        self.increase_count()
        self._count = type(self)._count

    @abstractmethod
    def __str__(self):
        """Return node content formatted for topology file"""

    def __repr__(self):
        """Default representation"""
        return f"{type(self).__name__}()"

    def _make_node_key(self) -> str:
        """Return string usable as node key"""
        return f"{self._node_key_name}_{self._count}"


class GenericNodeValue(NodeValue):
    """Generic fallback node value"""
    _node_key_name = "generic"

    def __init__(self, value,comment=None):
        super().__init__()
        self.value = value
        self.comment = comment

    def __str__(self):
        return self.value.__str__()

    def __repr__(self):
        return f"{type(self).__name__}(value={self.value})"


class Define(NodeValue):
    """#define or #undef directives"""
    _node_key_name = "define"

    def __init__(self, key, value, comment=None):
        super().__init__()
        self.key = key
        self.value = value
        self.comment = comment

    def __str__(self):
        p = f"  ; {self.comment}" if self.comment else ""
        if not isinstance(self.value, bool):
            return f"#define {self.key} {self.value}{p}"
        if self.value is True:
            return f"#define {self.key}{p}"
        if self.value is False:
            return f"#undef {self.key}{p}"

    def __repr__(self):
        return f"{type(self).__name__}(key={self.key}, value={self.value})"


class Condition(NodeValue):
    """#ifdef, #ifndef, #endif directives"""
    _node_key_name = "condition"

    def __init__(self, key, value, complement=None, comment=None):
        super().__init__()
        self.key = key
        self.value = value
        self.complement = complement
        self.comment = comment

    def __str__(self):
        p = f"  ; {self.comment}" if self.comment else ""
        if self.value is True:
            return f"#ifdef {self.key}{p}"
        if self.value is False:
            return f"#ifndef {self.key}{p}"
        if self.value is None:
            return "#endif"


class Section(NodeValue):
    """A regular section heading"""
    _node_key_name = "section"

    def __init__(self, title, comment=None):
        super().__init__()
        self.title = title
        self.comment = comment

    def __str__(self):
        if self.comment:
            return f"[{self.title}]  ; {self.comment}"
        return f"[{self.title}]"


class SpecializedSection(Section):
    category = 0
    allowed_occurrence = 0

    def __init__(self,comment=None):
        super().__init__(self._node_key_name,comment=comment)


class DefaultsSection(SpecializedSection):
    _node_key_name = "defaults"
    allowed_occurrence = 1


class AtomtypesSection(SpecializedSection):
    _node_key_name = "atomtypes"


class BondtypesSection(SpecializedSection):
    _node_key_name = "bondtypes"


class AngletypesSection(SpecializedSection):
    _node_key_name = "angletypes"


class PairtypesSection(SpecializedSection):
    _node_key_name = "pairtypes"


class DihedraltypesSection(SpecializedSection):
    _node_key_name = "dihedraltypes"


class ConstrainttypesSection(SpecializedSection):
    _node_key_name = "constrainttypes"


class NonbondedParamsSection(SpecializedSection):
    _node_key_name = "nonbonded_params"


class MoleculetypeSection(SpecializedSection):
    _node_key_name = "moleculetype"
    category = 1


class SystemSection(SpecializedSection):
    _node_key_name = "system"
    category = 2
    allowed_occurrence = 1


class MoleculesSection(SpecializedSection):
    _node_key_name = "molecules"
    category = 2
    allowed_occurrence = 1


class Subsection(Section):
    """A subsection heading"""
    _node_key_name = "subsection"

    def __init__(self, title, section=None, comment=None):
        super().__init__(title=title,comment=comment)
        self.section = section


class SpecializedSubsection(Subsection):
    category = 1
    allowed_occurrence = 0

    def __init__(self, section=None, comment=None):
        super().__init__(title=self._node_key_name,section=section,comment=comment)


class AtomsSubsection(SpecializedSubsection):
    _node_key_name = "atoms"


class BondsSubsection(SpecializedSubsection):
    _node_key_name = "bonds"


class PairsSubsection(SpecializedSubsection):
    _node_key_name = "pairs"


class PairsNBSubsection(SpecializedSubsection):
    _node_key_name = "pairs_nb"


class AnglesSubsection(SpecializedSubsection):
    _node_key_name = "angles"


class DihedralsSubsection(SpecializedSubsection):
    _node_key_name = "dihedrals"


class ExclusionsSubsection(SpecializedSubsection):
    _node_key_name = "exclusions"


class ConstraintsSubsection(SpecializedSubsection):
    _node_key_name = "constraints"


class SettlesSubsection(SpecializedSubsection):
    _node_key_name = "settles"


class VirtualSites2Subsection(SpecializedSubsection):
    _node_key_name = "virtual_sites2"


class VirtualSites3Subsection(SpecializedSubsection):
    _node_key_name = "virtual_sites3"


class VirtualSites4Subsection(SpecializedSubsection):
    _node_key_name = "virtual_sites4"


class VirtualSitesNSubsection(SpecializedSubsection):
    _node_key_name = "virtual_sitesn"


class PositionRestraintsSubsection(SpecializedSubsection):
    _node_key_name = "position_restraints"


class DistanceRestraintsSubsection(SpecializedSubsection):
    _node_key_name = "distance_restraints"


class DihedralRestraintsSubsection(SpecializedSubsection):
    _node_key_name = "dihedral_restraints"


class OrientationRestraintsSubsection(SpecializedSubsection):
    _node_key_name = "orientation_restraints"


class AngleRestraintsSubsection(SpecializedSubsection):
    _node_key_name = "angle_restraints"


class AngleRestraintsZSubsection(SpecializedSubsection):
    _node_key_name = "angle_restraints_z"


class Comment(NodeValue):
    """Standalone full-line comment"""
    _node_key_name = "comment"

    def __init__(self, comment: str):
        super().__init__()
        self.comment = comment
        self._char = ";"

    def __str__(self):
        if self._char is None:
            return f"{self.comment}"
        return f"{self._char} {self.comment.__str__()}"


class Include(NodeValue):
    """#include directive"""
    _node_key_name = "include"

    def __init__(self, include: str, comment=None):
        super().__init__()
        self.include = include
        self.comment = comment

    def __str__(self):
        p = f"  ; {self.comment}" if self.comment else ""
        return f"#include {self.include.__str__()}{p}"

    def __repr__(self):
        return f"{type(self).__name__}({self.include})"


class SectionEntry(NodeValue):
    """A section entry"""
    _node_key_name = "section_entry"
    _args = []

    def __init__(self, comment=None, **kwargs):
        super().__init__()
        for name, target_type, _ in self._args:
            value = kwargs.get(name)
            if value is None:
                setattr(self, name, None)
                continue

            if target_type is None:
                continue

            try:
                value = target_type(value)
            except ValueError:
                if not isinstance(value, str):
                    raise TypeError(
                        f"Argument {name} should be 'str' or convertible to {target_type.__name__}"
                    )
            setattr(self, name, value)

        self.comment = comment
        self._raw = None

    def __copy__(self):
        kwargs = {name:getattr(self,name) for name, *_ in self._args}
        copied = type(self)(**kwargs, comment=self.comment)
        copied._raw = self._raw
        return copied
    
    copy = __copy__     # alias

    @classmethod
    def create_raw(cls, line=None):
        node_value = cls()
        node_value._raw = line
        return node_value

    @classmethod
    def from_line(cls, *args, comment=None):
        kwargs = {kw[0]:v for kw,v in zip(cls._args, args)}
        entry = cls(comment=comment, **kwargs)
        return entry

    def __str__(self):
        if self._raw is not None:
            return self._raw
        result = ""
        for name, _, formatter in self._args:
            value = getattr(self, name)
            if value is None:
                continue
            if not isinstance(value, str) and isinstance(value, Iterable):
                for v in value:
                    result += f" {formatter(v)}"
                continue
            result += f" {formatter(value)}"
        if self.comment:
            result += f" ; {self.comment}"
        return result

    def __repr__(self):
        p = ", ".join(f"{name}={getattr(self,name)}" for name, *_ in self._args)
        return f"{type(self).__name__}({p})"


class PropertyInvoker:
    """Invokes descriptor protocol for instance attributes"""
    def __setattr__(self, attr, value):
        try:
            got = super().__getattribute__(attr)
            got.__set__(self, value)
        except AttributeError:
            super().__setattr__(attr, value)

    def __getattribute__(self, attr):
        got = super().__getattribute__(attr)
        try:
            return got.__get__(self, type(self))
        except AttributeError:
            return got


class P1TermEntry(SectionEntry, PropertyInvoker):
    _node_key_name = "p1_term_entry"
    _args = [
        ("i", int, _make_formatter(">6")),
        ("funct", int, _make_formatter(">6")),
        ("c", None, _make_formatter("12.6f")),
    ]

    def _getc(self, index=None):
        if index is None: index = 0
        def getter(self):
            return self.c[index]
        return getter

    def _setc(self, index=None):
        if index is None: index = 0
        def setter(self, value):
            self.c[index] = value
        return setter

    def _delc(self, index=None):
        if index is None: index = len(self.c) - 1
        def deleter(self):
            del self.c[index]
        return deleter

    def __init__(self,i=None,funct=None,c=None,comment=None,**kwargs):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_, **kwargs)
        if c is None: c = []
        self.c = [float(x) for x in c]
        for index in range(len(self.c)):
            # Add instance specific properties "c0", "c1", ...
            setattr(
                self, f"c{index}", property(
                    fget=self._getc(index),
                    fset=self._setc(index),
                    fdel=self._delc(index),
                    doc="Access term coefficients"
                )
            )

    @classmethod
    def from_line(cls, *args, comment=None):
        i, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        entry = cls(i=i, funct=funct, c=c, comment=comment)
        return entry


class P2TermEntry(P1TermEntry):
    _node_key_name = "p2_term_entry"
    _args = [
        ("i", int, _make_formatter(">6")),
        ("j", int, _make_formatter(">6")),
        ("funct", int, _make_formatter(">6")),
        ("c", None, _make_formatter("12.6f"))
    ]

    def __init__(self,i=None,j=None,funct=None,c=None,comment=None):
        super().__init__(i=i,j=j,funct=funct,c=c,comment=comment)

    @classmethod
    def from_line(cls, *args, comment=None):
        i, j, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        entry = cls(i=i,j=j,funct=funct,c=c,comment=comment)
        return entry


class P3TermEntry(P1TermEntry):
    _node_key_name = "p3_term_entry"
    _args = [
        ("i", int, _make_formatter(">6")),
        ("j", int, _make_formatter(">6")),
        ("k", int, _make_formatter(">6")),
        ("funct", int, _make_formatter(">6")),
        ("c", None, _make_formatter("12.6f"))
    ]

    def __init__(self,i=None,j=None,k=None,funct=None,c=None,comment=None):
        super().__init__(i=i,j=j,k=k,funct=funct,c=c,comment=comment)

    @classmethod
    def from_line(cls, *args, comment=None):
        i, j, k, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        entry = cls(i=i,j=j,k=k,funct=funct,c=c,comment=comment)
        return entry


class P4TermEntry(P1TermEntry):
    _node_key_name = "p4_term_entry"
    _args = [
        ("i", int, _make_formatter(">6")),
        ("j", int, _make_formatter(">6")),
        ("k", int, _make_formatter(">6")),
        ("l", int, _make_formatter(">6")),
        ("funct", int, _make_formatter(">6")),
        ("c", None, _make_formatter("12.6f"))
    ]

    def __init__(self,i=None,j=None,k=None,l=None,funct=None,c=None,comment=None):
        super().__init__(i=i,j=j,k=k,l=l,funct=funct,c=c,comment=comment)

    @classmethod
    def from_line(cls, *args, comment=None):
        i, j, k, l, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        entry = cls(i=i,j=j,k=k,l=l,funct=funct,c=c,comment=comment)
        return entry


class P6TermEntry(P1TermEntry):
    _node_key_name = "p4_term_entry"
    _args = [
        ("i", int, _make_formatter(">6")),
        ("j", int, _make_formatter(">6")),
        ("k", int, _make_formatter(">6")),
        ("l", int, _make_formatter(">6")),
        ("g", int, _make_formatter(">6")),
        ("h", int, _make_formatter(">6")),
        ("funct", int, _make_formatter(">6")),
        ("c", None, _make_formatter("12.6f"))
    ]

    def __init__(self,i=None,j=None,k=None,l=None,g=None,h=None,funct=None,c=None,comment=None):
        super().__init__(i=i,j=j,k=k,l=l,g=g,h=h,funct=funct,c=c,comment=comment)

    @classmethod
    def from_line(cls, *args, comment=None):
        i, j, k, l, g, h, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        entry = cls(i=i,j=j,k=k,l=l,g=g,h=h,funct=funct,c=c,comment=comment)
        return entry


class DefaultsEntry(SectionEntry):
    """Entry in defaults section"""
    _node_key_name = "defaults_entry"
    _args = [
        ("nbfunc", int, _make_formatter("<15")),
        ("comb_rule", int, _make_formatter("<15")),
        ("gen_pairs", str, _make_formatter("<15")),
        ("fudgeLJ", float, _make_formatter("<7")),
        ("fudgeQQ", float, _make_formatter("<7")),
        ("n", int, _make_formatter("<7"))
    ]

    def __init__(
            self,nbfunc=None,comb_rule=None,gen_pairs="no",fudgeLJ=None,fudgeQQ=None,
            n=None,comment=None
        ):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_)


class AtomtypesEntry(SectionEntry):
    _node_key_name = "atomtypes_entry"
    _args = [
        ("name", str, _make_formatter("<9")),
        ("bond_type", str, _make_formatter("<4")),
        ("at_num", int, _make_formatter("<3")),
        ("mass", float, _make_formatter("<8")),
        ("charge", float, _make_formatter("<6")),
        ("ptype", str, _make_formatter("<1")),
        ("sigma", float, _make_formatter("1.5e")),
        ("epsilon", float, _make_formatter("1.5e"))
    ]

    def __init__(
            self, name=None, bond_type=None, at_num=None, mass=None,
            charge=None, ptype=None, sigma=None, epsilon=None, comment=None
        ):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_)

    @classmethod
    def from_line(cls, *args, comment):
        if len(args) == 7:          # special case
            arg_names = cls._args[0:1] + cls._args[2:]
        else:
            arg_names = cls._args
        kwargs = {kw[0]:v for kw,v in zip(arg_names, args)}
        entry = cls(comment=comment,**kwargs)
        return entry


class BondtypesEntry(P2TermEntry):
    _node_key_name = "bondtypes_entry"


class AngletypesEntry(P3TermEntry):
    _node_key_name = "angletypes_entry"


class PairtypesEntry(P2TermEntry):
    _node_key_name = "pairtypes_entry"


class DihedraltypesEntry(P4TermEntry):
    _node_key_name = "dihedraltypes_entry"


class ConstrainttypesEntry(P2TermEntry):
    _node_key_name = "constrainttypes_entry"


class NonbondedParamsEntry(P2TermEntry):
    _node_key_name = "nonbonded_params_entry"


class MoleculetypeEntry(SectionEntry):
    _node_key_name = "moleculetype_entry"
    _args = [
        ("molecule", str, _make_formatter("")),
        ("nrexcl", int, _make_formatter(">6"))
    ]

    def __init__(self,molecule=None,nrexcl=None,comment=None):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_)


class SystemEntry(SectionEntry):
    _node_key_name = "system_entry"
    _args = [
        ("name", str, _make_formatter(""))
    ]

    def __init__(self,name=None,comment=None):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_)

    @classmethod
    def from_line(cls, *args, comment):
        entry = cls(name=" ".join(args),comment=comment)
        return entry


class MoleculesEntry(SectionEntry):
    _node_key_name = "molecules_entry"
    _args = [
        ("molecule", str, _make_formatter("")),
        ("number", int, _make_formatter(">6"))
    ]

    def __init__(self,molecule=None,number=None,comment=None):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_)


class AtomsEntry(SectionEntry):
    _node_key_name = "atoms_entry"
    _args = [
        ("nr", int, _make_formatter("<6")),
        ("type", str, _make_formatter("<6")),
        ("resnr", int, _make_formatter("<6")),
        ("residue", str, _make_formatter("<6")),
        ("atom", str, _make_formatter("<6")),
        ("cgnr", int, _make_formatter("<6")),
        ("charge", float, _make_formatter(">10.6f")),
        ("mass", float, _make_formatter("8.3f")),
        ("typeB", float, _make_formatter("<6")),
        ("chargeB", float, _make_formatter(">10.6f")),
        ("massB", float, _make_formatter("8.3f"))
    ]

    def __init__(
            self, nr=None, type=None, resnr=None, residue=None,
            atom=None, cgnr=None, charge=None, mass=None,
            typeB=None, chargeB=None, massB=None, comment=None
        ):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_)


class BondsEntry(P2TermEntry):
    _node_key_name = "bonds_entry"


class PairsEntry(P2TermEntry):
    _node_key_name = "pairs_entry"


class PairsNBEntry(P2TermEntry):
    _node_key_name = "pairs_nb_entry"


class AnglesEntry(P3TermEntry):
    _node_key_name = "angles_entry"


class DihedralsEntry(P4TermEntry):
    _node_key_name = "dihedrals_entry"


class ExclusionsEntry(SectionEntry):
    _node_key_name = "exclusions_entry"
    _args = [
        ("indices", None, _make_formatter(">6"))
    ]

    def __init__(self,indices=None,comment=None):
        locals_ = _trim_locals(locals())
        super().__init__(**locals_)
        if indices is None: indices = []
        self.indices = [int(i) for i in indices]

    @classmethod
    def from_line(cls, *args, comment):
        entry = cls(indices=args,comment=comment)
        return entry


class ConstraintsEntry(P2TermEntry):
    _node_key_name = "constraints_entry"


class SettlesEntry(P1TermEntry):
    _node_key_name = "settles_entry"


class VirtualSites1Entry(P1TermEntry):
    _node_key_name = "virtual_sites1_entry"
    _args = [
        ("i", int, _make_formatter(">6")),
        ("f", None, _make_formatter(">6")),
        ("funct", int, _make_formatter(">6")),
        ("c", None, _make_formatter("12.6f"))
    ]

    def __init__(self,i=None,f=None,funct=None,c=None,comment=None):
        super().__init__(i=i,funct=funct,c=c,comment=comment)
        if f is None: f = []
        self.f = [int(x) for x in f]

    @classmethod
    def from_line(cls, *args, comment=None):
        i, f1, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        f = [f1]
        entry = cls(i=i,f=f,funct=funct,c=c,comment=comment)
        return entry


class VirtualSites2Entry(VirtualSites1Entry):
    _node_key_name = "virtual_sites2_entry"

    @classmethod
    def from_line(cls, *args, comment=None):
        i, f1, f2, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        f = [f1, f2]
        entry = cls(i=i,f=f,funct=funct,c=c,comment=comment)
        return entry


class VirtualSites3Entry(VirtualSites1Entry):
    _node_key_name = "virtual_sites3_entry"

    @classmethod
    def from_line(cls, *args, comment=None):
        i, f1, f2, f3, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        f = [f1, f2, f3]
        entry = cls(i=i,f=f,funct=funct,c=c,comment=comment)
        return entry


class VirtualSites4Entry(VirtualSites1Entry):
    _node_key_name = "virtual_sites4_entry"

    @classmethod
    def from_line(cls, *args, comment=None):
        i, f1, f2, f3, f4, *rest = args
        if rest:
            funct, *c = rest
        else:
            funct = c = None
        f = [f1, f2, f3, f4]
        entry = cls(i=i,f=f,funct=funct,c=c,comment=comment)
        return entry


class VirtualSitesNEntry(VirtualSites1Entry):
    _node_key_name = "virtual_sitesn_entry"
    _args = [
        ("i", int, _make_formatter(">6")),
        ("funct", int, _make_formatter(">6")),
        ("f", None, _make_formatter(">6")),
    ]

    def __init__(self,i=None,funct=None,f=None,comment=None):
        super().__init__(i=i,funct=funct,f=f,comment=comment)
        delattr(self, "c")

    @classmethod
    def from_line(cls, *args, comment=None):
        i, *rest = args
        if rest:
            funct, *f = rest
        else:
            funct = f = None
        entry = cls(i=i,f=f,funct=funct,comment=comment,)
        return entry


class PositionRestraintsEntry(P1TermEntry):
    _node_key_name = "position_restraints_entry"


class DistanceRestraintsEntry(P2TermEntry):
    _node_key_name = "distance_restraints_entry"


class DihedralRestraintsEntry(P4TermEntry):
    _node_key_name = "dihedral_restraints_entry"


class OrientationRestraintsEntry(P2TermEntry):
    _node_key_name = "orientations_restraints_entry"


class AngleRestraintsEntry(P4TermEntry):
    _node_key_name = "angle_restraints_entry"


class AngleRestraintsZEntry(P2TermEntry):
    _node_key_name = "angle_restraints_z_entry"


DEFAULT_NODES = {
    "generic": GenericNodeValue,
    "comment": Comment,
    "define": Define,
    "include": Include,
    "condition": Condition,
    "section": Section,
    "subsection": Subsection,
    "entry": SectionEntry,
    "virtual_sites1_entry": VirtualSites1Entry,
    "p1term_entry": P1TermEntry,
    "p2term_entry": P2TermEntry,
    "p3term_entry": P3TermEntry,
    "p4term_entry": P4TermEntry,

    "atomtypes": AtomtypesSection,              "atomtypes_entry": AtomtypesEntry,
    "bondtypes": BondtypesSection,              "bondtypes_entry": BondtypesEntry,
    "angletypes": AngletypesSection,            "angletypes_entry": AngletypesEntry,
    "pairtypes": PairtypesSection,              "pairtypes_entry": PairtypesEntry,
    "dihedraltypes": DihedraltypesSection,      "dihedraltypes_entry": DihedraltypesEntry,
    "constrainttypes": ConstrainttypesSection,  "constrainttypes_entry": ConstrainttypesEntry,
    "moleculetype": MoleculetypeSection,        "moleculetype_entry": MoleculetypeEntry,
    "system": SystemSection,                    "system_entry": SystemEntry,
    "molecules": MoleculesSection,              "molecules_entry": MoleculesEntry,
    "defaults": DefaultsSection,                "defaults_entry": DefaultsEntry,
    "nonbonded_params": NonbondedParamsSection, "nonbonded_params_entry": NonbondedParamsEntry,
    "atoms": AtomsSubsection,                   "atoms_entry": AtomsEntry,
    "bonds": BondsSubsection,                   "bonds_entry": BondsEntry,
    "pairs": PairsSubsection,                   "pairs_entry": PairsEntry,
    "pairs_nb": PairsNBSubsection,              "pairs_nb_entry": PairsNBEntry,
    "angles": AnglesSubsection,                 "angles_entry": AnglesEntry,
    "dihedrals": DihedralsSubsection,           "dihedrals_entry": DihedralsEntry,
    "exclusions": ExclusionsSubsection,         "exclusions_entry": ExclusionsEntry,
    "constraints": ConstraintsSubsection,       "constraints_entry": ConstraintsEntry,
    "settles": SettlesSubsection,               "settles_entry": SettlesEntry,
    "virtual_sites2": VirtualSites2Subsection,  "virtual_sites2_entry": VirtualSites2Entry,
    "virtual_sites3": VirtualSites3Subsection,  "virtual_sites3_entry": VirtualSites3Entry,
    "virtual_sites4": VirtualSites4Subsection,  "virtual_sites4_entry": VirtualSites4Entry,
    "virtual_sitesn": VirtualSitesNSubsection,  "virtual_sitesn_entry": VirtualSitesNEntry,
    "position_restraints": PositionRestraintsSubsection, "position_restraints_entry": PositionRestraintsEntry,
    "distance_restraints": DistanceRestraintsSubsection, "distance_restraints_entry": DistanceRestraintsEntry,
    "dihedral_restraints": DihedralRestraintsSubsection, "dihedral_restraints_entry": DihedralRestraintsEntry,
    "angle_restraints": AngleRestraintsSubsection,       "angle_restraints_entry": AngleRestraintsEntry,
    "angle_restraints_z": AngleRestraintsZSubsection,    "angle_restraints_z_entry": AngleRestraintsZEntry,
    "orientation_restraints": OrientationRestraintsSubsection,
    "orientation_restraints_entry": OrientationRestraintsEntry,
}

SUPPORTED_DIRECTIVES = sorted([
    k for k in DEFAULT_NODES.keys() if not (k.endswith('_entry') or k in [
        'generic','comment','condition','section','subsection','entry',
    ])
])


class Node:
    __slots__ = ["prev", "next", "key", "value", '__weakref__']

    def __init__(self):
        self.prev = None
        self.next = None
        self.key = None
        self.value = None

    def __repr__(self):
        attr_str = f"(key={self.key}, value={self.value})"
        return f"{type(self).__name__}{attr_str}"

    def connect(self, other, forward=True):
        """Link another node in forward/backward direction"""
        if forward is True:
            self.next = other
            other.prev = weakref.proxy(self)
        else:
            self.prev = weakref.proxy(other)
            other.next = self


class GromacsTop:
    __node_value_types = DEFAULT_NODES

    def __init__(self):
        self._nodes = dict()
        self._hardroot = Node()     # the first `node` will always be empty `Node`
        self._root = root = weakref.proxy(self._hardroot)
        root.prev = root.next = root

    def __str__(self):
        result = ""
        for node in self:
            if isinstance(node.value, self.__node_value_types["section"]):
                result += "\n"
            result += f"{node.value}\n"
        return result

    def __repr__(self):
        return f"{type(self).__name__}()"

    def __iter__(self):
        root = current = self._root
        while current.next is not root:     # `not None` to avoid cyclic link
            current = current.next
            yield current

    def __reversed__(self):
        root = current = self._root
        while current.prev is not root:
            current = current.prev
            yield current

    def __len__(self):
        """Return number of linked nodes"""
        length = 0
        for length, _ in enumerate(self, 1):
            pass
        return length

    def __getitem__(self, query):
        if isinstance(query, str):
            return self._nodes.__getitem__(query)
        if isinstance(query, slice):
            return itertools.islice(self, query.start, query.stop, query.step)
        if isinstance(query, int):
            if query < 0:
                query = (query * -1) - 1
                iterable = reversed(self)
            else:
                iterable = iter(self)
            try:
                return next(itertools.islice(iterable, query, query + 1))
            except StopIteration:
                raise IndexError("index out of range")
        raise ValueError(
            f"items can only be queried by 'str' (node key), "
            f"or 'int' or 'slice' (node index), not {type(query).__name__}"
        )

    def __contains__(self, key):
        if key in self._nodes:
            return True
        return False

    def _check_key_and_add_new_node(self, key):
        if key in self._nodes:
            raise KeyError(f"node {key} already exist")
        self._nodes[key] = node = Node()
        return node

    @classmethod
    def make_nvtype(cls, name, *args, **kwargs):
        """Retrieve node type by name, initialize, and make key"""
        nvtype = cls.__node_value_types[name]
        node_value = nvtype(*args, **kwargs)
        node_key = node_value._make_node_key()
        return node_key, node_value

    @classmethod
    def select_nvtype(cls, name):
        """Retrieve node type by name"""
        return cls.__node_value_types[name]

    def add(self, key, value) -> None:
        node = self._check_key_and_add_new_node(key)
        root = self._root
        last = root.prev
        node.prev, node.next, node.key, node.value = last, root, key, value  # direction linked list
        last.next = node
        root.prev = weakref.proxy(node)

    def pop(self, key):
        node = self._nodes.pop(key)
        node.prev.next = node.next
        node.next.prev = node.prev
        return node

    def discard(self, key) -> None:
        try:
            node = self._nodes.pop(key)
        except KeyError:
            pass
        else:
            node.prev.next = node.next
            node.next.prev = node.prev

    def replace(self, key, value):
        """Replace node with specific key while retaining key"""
        old_node = self._nodes.pop(key)
        self._nodes[key] = new_node = Node()
        new_node.key, new_node.value = key, value
        new_node.prev, new_node.next = old_node.prev, old_node.next
        new_node.prev.next = new_node
        new_node.next.prev = weakref.proxy(new_node)

    def index(self, key, start=None, stop=None):
        """Return index of node

        Args:
            key: Node key
            start: Ignore nodes with lower index
            stop: Ignore nodes with index greater or equal

        Raises:
            ValueError if the key is not present
        """
        if not start: start = -1
        for i, node in enumerate(self):
            if i < start:
                continue
            if stop and i >= stop:
                break
            if node.key == key:
                return i
        raise ValueError(f"node {key} is not in list")

    def insert(self, index, key, value):
        """Insert node before index"""
        node = self._check_key_and_add_new_node(key)
        prev = root = self._root
        next_node = root.next
        for i, next_node in enumerate(self):
            if i == index:
                prev = next_node.prev
                break
        else:
            prev = root.prev
            next_node = prev.next
        node.prev, node.next = prev, next_node
        node.key, node.value = key, value
        prev.next = node
        next_node.prev = weakref.proxy(node)

    def relative_insert(self, node, key, value, forward=True):
        """Insert node after/before other node"""
        new_node = self._check_key_and_add_new_node(key)
        if forward is True:
            new_node.prev, new_node.next = node.next.prev, node.next
            new_node.next.prev = weakref.proxy(new_node)
            node.next = new_node
        else:
            new_node.prev, new_node.next = node.prev, node
            new_node.prev.next = new_node
            node.prev = weakref.proxy(new_node)
        new_node.key, new_node.value = key, value

    def block_insert(self, node, block_start, block_end, forward=True):
        """Insert block of linked nodes before/after node"""
        if forward is True:
            block_start.prev = node.next.prev
            block_end.next = node.next
            block_end.next.prev = weakref.proxy(block_end)
            node.next = block_start
        else:
            block_start.prev = node.prev
            block_end.next = node
            node.prev.next = block_start
            node.prev = weakref.proxy(block_end)
        current = block_start
        while current is not block_end.next:
            if current.key in self._nodes:
                raise KeyError(f"node {current.key} does already exist")
            self._nodes[current.key] = current
            current = current.next

    def block_discard(self, block_start, block_end):
        block_start.prev.next = block_end.next
        block_end.next.prev = block_start.prev
        block_start.prev = None
        block_end.next = None
        current = block_start
        while True:
            _ = self._nodes[current.key]
            if current is block_end:
                break
            current = current.next

    def get_next_node_with_nvtype(self,start=None,stop=None,nvtype=None,exclude=None,forward=True):
        """Search topology for another node

        Args:
            start(Node): If `None`, an `nvtype` must exist and search will start from the beginning.
            stop(Node) : If `None`, search until end.
            nvtype     : Type of node value to search for. If `None`, search for same type as start.
            exclude    : Exclude node types from search.
            forward    : If `True`, search topology forwards.  If `False`, search backwards.
        """
        if start is None:
            if nvtype is None:
                raise ValueError("If start=None, a node type must be specified")
            start = self._root

        if stop is None: stop = self._root
        if nvtype is None: nvtype = type(start.value)
        if exclude is None: exclude = ()

        goto = "next" if forward else "prev"
        node = getattr(start, goto)
        while node.prev.next is not stop:
            if not isinstance(node.value, nvtype) or isinstance(node.value, exclude):
                node = getattr(node, goto)
            else:
                return node
        raise LookupError(f"Node of type {nvtype} not found")  # this should not happen

    @property
    def includes_resolved(self):
        for node in self:
            if isinstance(node.value, self.__node_value_types["include"]):
                return False
        return True

    @property
    def conditions_resolved(self):
        for node in self:
            if isinstance(node.value, self.__node_value_types["condition"]):
                return False
        return True

    def find_complement(self, node):
        """Find complementary Condition node"""
        root = self._root
        if node.value.value is None:
            goto = "prev"
        else:
            goto = "next"
        current = getattr(node, goto)
        while current is not root:
            if isinstance(current.value, self.__node_value_types["condition"]):
                if current.value.key == node.value.key:
                    return current
            current = getattr(current, goto)

    def sanitize(self):
        """remove empty section, and return its values"""
        keys = []
        stop = self._root
        node = self._root.next
        while node is not stop:
            if isinstance(node.value,Subsection):
                key = node.key
                while node is not stop:
                    node = node.next
                    if isinstance(node.value,Comment):
                        continue
                    break
                if not isinstance(node.value,SectionEntry):
                    keys.append(key)
            else:
                node = node.next
        for k in keys: self.pop(k)
        return keys


class GromacsTopParser:
    """Read and write GROMACS topology files
    
    Args:
        parse_includes(bool)    : whether parse `#include`
        custom_paths(list)      : where to find file for `#include`
        include_blacklist(list) : ignore inside files for `#include`

        definitions(dict)       : conditions for `#ifdef`
        resolve_conditions(bool): whether resolve `#ifdef`, if defined, parse; else, skip
    """
    __top_type = GromacsTop

    def __init__(
            self, parse_includes=True,custom_paths=None,include_blacklist=None,
            definitions=None,resolve_conditions=True,
        ):
        self.parse_includes = True if parse_includes is True else False
        if custom_paths:
            if isinstance(custom_paths,str):
                self.custom_paths = [custom_paths, ]
            elif isinstance(custom_paths,(list,tuple)):
                self.custom_paths = [i for i in custom_paths if isinstance(i,str)]
            else:
                print(f'Warning: wrong defined: custom_paths: {custom_paths} => ignoring')
                self.custom_paths = []
        else:
            self.custom_paths = []
        if include_blacklist:
            if isinstance(include_blacklist,str):
                self.include_blacklist = [include_blacklist, ]
            elif isinstance(include_blacklist,(list,tuple)):
                self.include_blacklist = [i for i in include_blacklist if isinstance(i,str)]
            else:
                print(f'Warning: wrong defined: include_blacklist: {include_blacklist} => ignoring')
                self.include_blacklist = []
        else:
            self.include_blacklist = []

        self.resolve_conditions = resolve_conditions
        if definitions and isinstance(definitions,dict):
            self.definitions.update(definitions)
        else:
            if definitions is not None:
                print(f'Warning: wrong defined: definitions(dict): {definitions} => ignoring')
            self.definitions = {}

    def _preprocess_includes(
            self, file,     # file can be literally anything that will be correctly parsed
            parse_includes=None,custom_paths=None,include_blacklist=None,
            processed_paths=None,parsed_full_paths=None,
        ):
        """Pre-process topology file line by line and resolve '#include'"""
        if hasattr(file,'name'):
            pass
        elif isinstance(file,str) and os.path.isfile(file):
            file = open(file,'rt')      # GC will do `file.close`
        elif isinstance(file,io.StringIO):
            if not hasattr(file,'name'): file.name = 'io.StringIO'
        else:
            return '',None

        if not processed_paths: processed_paths = []
        if file.name in processed_paths:
            print(f'  => Warning: double include: {file.name}  => ignoring')
            return '',None  # to aviod infinite loop or double parsing
        processed_paths.append(file.name)

        if not parse_includes:
            previous = ''
            for line in file:
                if line.endswith(' \\'):
                    # Resolve multi-line statement
                    previous = f"{previous}{line[:-1]}"  # tail space will always be there
                    continue
                line = f"{previous}{line}"
                previous = ''
                yield self.split_comment(line)
            return '',None

        if parsed_full_paths is None:
            if custom_paths:
                parsed_full_paths = [('(custom)',pathlib.Path(p)) for p in custom_paths]
            else:
                parsed_full_paths = []
                try:
                    file_path = pathlib.Path(file.name).parent.absolute()
                except AttributeError:
                    pass
                else:
                    parsed_full_paths.append(('(local)',file_path))
                p = self.get_gmx_datadir()
                if p:
                    parsed_full_paths.append(('(shared)',p))

        previous = ''
        for line in file:
            if line.endswith(' \\'):
                # Resolve multi-line statement
                previous = f"{previous}{line[:-1]}"  # tail space will always be there
                continue
            line = f"{previous}{line}"
            previous = ''

            line,comment = self.split_comment(line)
            if not line.startswith('#include'):
                yield line,comment
                continue

            found = False
            include_file = line.split()[1].strip('"').strip()
            for label,include_dir in parsed_full_paths:
                include_path = include_dir / include_file
                if not include_path.is_file():
                    continue

                if include_blacklist:
                    g = str(include_path)
                    for t in include_blacklist:
                        if t in g:
                            found = True
                            print(f"Not including {include_path} (blacklist)")
                            break
                    if found:
                        break

                print(f"Including {include_path} {label}")
                with open(include_path) as f:
                    yield from self._preprocess_includes(
                        f,include_blacklist=include_blacklist,
                        processed_paths=processed_paths, parsed_full_paths=parsed_full_paths
                    )
                found = True
                break

            if not found:
                print(f"Could not find {include_file}")
                yield line,comment

    def read(self, file) -> GromacsTop:
        file = self._preprocess_includes(
            file, parse_includes=self.parse_includes,custom_paths=self.custom_paths,
            include_blacklist=self.include_blacklist
        )

        top = self.__top_type()
        for node_value_type in top._GromacsTop__node_value_types.values():
            node_value_type.reset_count()

        active_section = None
        active_category = 0
        active_conditions = OrderedDict()
        active_definitions = {}
        active_definitions.update(self.definitions)
        for line,comment in file:
            # line = line.strip()       # already happened
            if not line:
                if comment:
                    node_key, node_value = top.make_nvtype("comment", comment)
                    top.add(node_key, node_value)
                continue

            if line.startswith('#define'):
                line = line.lstrip("#define").split(maxsplit=1)
                if len(line) == 1:
                    node_key, node_value = top.make_nvtype("define", line[0], True, comment=comment)
                    top.add(node_key, node_value)
                    active_definitions[line[0]] = True
                else:
                    node_key, node_value = top.make_nvtype("define", line[0], line[1], comment=comment)
                    top.add(node_key, node_value)
                    active_definitions[line[0]] = line[1]
                continue

            if line.startswith('#undef'):
                line = line.lstrip("#undef").lstrip()
                node_key, node_value = top.make_nvtype("define", line, False, comment=comment)
                top.add(node_key, node_value)
                _ = active_definitions.pop(line)
                continue

            if line.startswith('#ifdef'):
                line = line.lstrip('#ifdef').lstrip()
                active_conditions[line] = True
                if not self.resolve_conditions:
                    node_key, node_value = top.make_nvtype("condition", line, True, comment=comment)
                    top.add(node_key, node_value)
                continue

            if line.startswith('#ifndef'):
                line = line.lstrip('#ifndef').lstrip()
                active_conditions[line] = False
                if not self.resolve_conditions:
                    node_key, node_value = top.make_nvtype("condition", line, False, comment=comment)
                    top.add(node_key, node_value)
                continue

            if line.startswith('#else'):
                last_condition, last_value = next(reversed(active_conditions.items()))
                rev = not last_value
                active_conditions[last_condition] = rev

                if not self.resolve_conditions:
                    node_key, node_value = top.make_nvtype("condition", last_condition, None, comment=comment)
                    top.add(node_key, node_value)

                    node_key, node_value = top.make_nvtype("condition", last_condition, rev, comment=comment)
                    top.add(node_key, node_value)
                    continue

            if line.startswith('#endif'):
                last_condition, _ = active_conditions.popitem(last=True)
                if not self.resolve_conditions:
                    node_key, node_value = top.make_nvtype("condition", last_condition, None, comment=comment)
                    top.add(node_key, node_value)

                    node = top[-1]
                    complement = top.find_complement(node)
                    if complement is not None:
                        node.value.complement = self.ensure_proxy(complement)
                        complement.value.complement = self.ensure_proxy(node)
                continue

            if self.resolve_conditions:
                skip = False
                for condition, required_value in active_conditions.items():
                    defined_value = active_definitions.get(condition, False)
                    # `required_value` will always be `True`
                    # thus, no need check `defined_value == required_value`
                    if not defined_value:
                        skip = True
                        break
                if skip:
                    continue

            if line.startswith("#include"):
                include = line.strip("#include").lstrip()
                node_key, node_value = top.make_nvtype("include", include, comment=comment)
                top.add(node_key, node_value)
                continue

            if line.startswith('['):
                _new_section = line.strip(' []').casefold()
                nvtype = top._GromacsTop__node_value_types.get(_new_section, None)
                if nvtype is None:
                    # Should not happen for compliant topologies
                    print(f"Unknown section {_new_section}")
                    node_key, node_value = top.make_nvtype("section", _new_section, comment=comment)
                    top.add(node_key, node_value)
                    active_section = node_value
                    continue

                if (nvtype.category < active_category):
                    print(f"Inconsistent section {_new_section}")
                else:
                    active_category = nvtype.category

                issubsection = issubclass(nvtype, top._GromacsTop__node_value_types["subsection"])
                if issubsection:
                    node_value = nvtype(section=weakref.proxy(active_section),comment=comment)
                    active_section = node_value
                else:
                    node_value = nvtype(comment=comment)
                    active_section = node_value

                node_key = node_value._make_node_key()
                top.add(node_key, node_value)
                continue

            if active_section is None:
                node_key, node_value = top.make_nvtype("comment", line, comment=comment)
                node_value._char = None
                top.add(node_key, node_value)
                continue

            expected_entry = f"{active_section._node_key_name}_entry"
            nvtype = top._GromacsTop__node_value_types.get(expected_entry, None)
            if nvtype:
                args = line.split()
                try:
                    node_value = nvtype.from_line(*args, comment=comment)
                except TypeError:
                    # Should not happen if entry type can deal with line
                    if comment: line += f" ; {comment}"
                    node_value = nvtype.create_raw(f"{line}")
                finally:
                    node_key = node_value._make_node_key()
                top.add(node_key, node_value)
                continue

            # Absolute fallback
            node_key, node_value = top.make_nvtype("generic", line)
            top.add(node_key, node_value)

        return top

    def get_gmx_datadir(self):
        """actually `stderr` will be used"""
        calls = [['gmx','-h'], ['gmx_mpi','-h'], ['gmx_cpu','-h']]
        for c in calls:
            try:
                p = subprocess.run(c,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf8')
            except FileNotFoundError:
                return
            else:
                if p.returncode == 0:
                    break
                else:
                    return
        _feedback = p.stderr.split('\n') + p.stdout.split('\n')
        gmx_shared = None
        for line in _feedback:
            # if line.startswith('Executable:'):
            #     gmx_exe = pathlib.Path(line.split()[-1])
            if line.startswith('Data prefix:'):
                gmx_shared = pathlib.Path(line.split()[-1]) / 'share/gromacs/top'
                break
        return gmx_shared

    def split_comment(self,line):
        n = line.find(';')
        if n != -1:
            return line[:n].strip(), line[n+1:].strip()
        return line.strip(), None

    def ensure_proxy(self,obj):
        """Return a proxy of an object avoiding proxy of proxy"""
        if not isinstance(obj, (weakref.ProxyType, weakref.CallableProxyType)):
            return weakref.proxy(obj)
        return obj


class GromacsTopSelection:
    def __init__(self, top=None, *args, **kws):
        self.top = top if isinstance(top,GromacsTop) else None
        self._directive_keys = []
        self._molecules_keys = []
        self._molecules_entry_keys = []

    def get_all_parsed_directive_keys(self):
        if not self.top: return self._directive_keys
        self._directive_keys = []
        section_nvtype = GromacsTop.select_nvtype("section")
        for node in self.top:
            if isinstance(node.value, section_nvtype):
                self._directive_keys.append(node.key)
        return self._directive_keys

    def get_directive_node_from_key(self,key):
        if not self._directive_keys:
            self.get_all_parsed_directive_keys()
        if not self._directive_keys:
            print('Warning: emptyp directive keys')
            return
        if key not in self._directive_keys:
            print(f'Warning: key not exist: {key}')
            return
        for node in self.top:
            if node.key == key:
                return node
            node = node.next

    def get_directive_nodes_from_directive(self,directive):
        if not self.top: return
        if directive not in SUPPORTED_DIRECTIVES:
            raise NotImplementedError('directive not implemented')
        nvtype = GromacsTop.select_nvtype(directive)
        results = []
        for node in self.top:
            if isinstance(node.value, nvtype):
                results.append(node)
        return results

    def get_all_lines_starting_comments(self):
        """get all comment lines like: `; comment at the very beginning`"""
        results = self.get_all_lines_contain_comments()
        return [l for l in results if l.startswith(';')]
    

    def get_all_lines_contain_comments(self):
        """get all comment lines that have comments`"""
        if not self.top: return
        results = []
        for node in self.top:
            if node.value.comment:
                results.append(node.value.__str__())
        return results

    def get_all_lines_inline_comments(self):
        """get all comment lines like: `[system] ; comment after`"""
        results = self.get_all_lines_contain_comments()
        return [l for l in results if not l.startswith(';')]

    def get_all_molecules_type_keys(self):
        if not self.top: return
        self._molecules_keys = []
        for node in self.top:
            if isinstance(node.value, MoleculetypeEntry):
                self._molecules_keys.append(node.value.molecule)
        return self._molecules_keys

    def get_all_molecules_entry_keys(self):
        if not self.top: return
        if self._molecules_entry_keys: return
        self._molecules_entry_keys = []
        for node in self.top:
            if isinstance(node.value,Subsection):
                self._molecules_entry_keys.append(node.key)
        return self._molecules_entry_keys

    def get_all_atoms_and_bonds_sections(self):
        if not self.top: return
        sections = {}
        stop = self.top._root
        node = self.top._root.next
        while node is not stop:
            if isinstance(node.value,Subsection):
                key = node.key
                while node is not stop:
                    node = node.next
                    if isinstance(node.value,Comment):
                        continue
                    if isinstance(node.value,SectionEntry):
                        if 'atoms' in key:
                            vl = sections.setdefault(key, [])
                            vl.append([
                                node.value.nr, node.value.type, node.value.resnr, node.value.residue,
                                node.value.atom, node.value.cgnr, node.value.charge, node.value.mass,
                                node.value.typeB, node.value.chargeB, node.value.massB,
                            ])
                        elif 'bonds' in key:
                            vl = sections.setdefault(key, [])
                            vl.append([
                                node.value.i, node.value.j, node.value.funct, *node.value.c
                            ])
                    else:
                        break
            else:
                node = node.next
        return sections


class GromacsTopOps(GromacsTopSelection):
    """
    Args:
        sections(dict):
            {'directive':List[pars], 'comment':str, 'gro':List[pars], 'xyz':List[float] }

        >> gro: List[[resnum, resname, atomname, x, y, z, vx, vy, vz]] : (x:8.3f, vx:8.4f)
        >> xyz: List[[atomname, x, y, z]]
    """
    def __init__(self,topfile=None,sections=None,*args,**kws):
        top = kws.get('top',None)
        if top:
            if topfile:
                print('Warning: double definition: top & topfile: top will be used instead')
        else:
            if topfile:
                parser = GromacsTopParser()
                kws['top'] = parser.read(topfile)
            else:
                print('Warning: not defined: top')
        super().__init__(*args,**kws)
        self.kws = kws
        self.sections = sections if sections else {}
        if self.top:
            sec = self.get_all_atoms_and_bonds_sections()
            if sec:
                self.sections.update(sec)
                self.sections['key_id'] = list(sec.keys())
        pt = Chem.GetPeriodicTable()
        self._periodic_table = [pt.GetElementSymbol(i) for i in range(101)]

    def get_mol(self):
        if 'gro' in self.sections:
            mol = self.get_mol_from_gro()
        elif 'xyz' in self.sections:
            mol = self.get_mol_from_xyz()
        else:
            mol = self.get_mol_from_top()
        return mol

    def get_mol_from_top(self,key_id=None):
        """key_id: integer"""
        print('Warning: TODO: not yet done')
        return 
        if key_id:
            aid = 'atoms_' + str(key_id)
            bid = 'bonds_' + str(key_id)
            if aid in self.sections and bid in self.sections:
                secatoms = self.sections[aid]
                secbonds = self.sections[bid]
            else:
                print(f'Warning: not exist key_id: {key_id}')
                print(f'  --> found key_id: {self.sections["key_id"]}')
                return
        elif 'atoms' in self.sections and 'bonds' in self.sections:
            secatoms = self.sections['atoms']
            secbonds = self.sections['bonds']
        else:
            return
        from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
        mol = Chem.RWMol()
        for g in secatoms:
            a = self.guess_elemtype(g[4])
            mol.AddAtom(Chem.Atom(a))
        bonds = []
        for g in secbonds:
            mol.AddBond(g[0]-1,g[1]-1)
            bonds.append((g[0]-1,g[1]-1,g[3]*10))   # unit nanometer to angstrom
        mol = mol.GetMol()
        if mol.GetNumAtoms() <= 10:
            n = 20
        elif mol.GetNumAtoms() <= 30:
            n = 50
        else:
            n = 100
        EmbedMultipleConfs(mol,numConfs=n,maxAttempts=500)
        rmsd = []
        for conf in mol.GetConformers():
            tot = 0.0
            for p in bonds:
                d = conf.GetAtomPosition(p[0]).Distance(conf.GetAtomPosition(p[1]))
                tot += (d-p[2]) ** 2
            rmsd.append(tot)
        t = rmsd.index(min(rmsd))
        sel = mol.GetConformer(t)
        new = Chem.RWMol(mol)
        new.RemoveAllConformers()
        new.AddConformer(sel)
        mol = new.GetMol()
        if 'moleculetype' in self.sections:
            name = self.sections['moleculetype'][0][0]
            mol.SetProp('_Name',name)
        mol.SetProp('Title',self.kws['title']) if 'title' in self.kws else mol.SetProp('Title','GMXTop')
        for i,a in enumerate(mol.GetAtoms()):
            a.SetProp('name',secatoms[i][4])
            a.SetProp('residue',secatoms[i][3])
        return mol

    def get_mol_from_gro(self,section=None):
        if not section:
            if not 'gro' in self.sections: return
            section = self.sections['gro']
        pdb = ''
        for i,g in enumerate(section):
            r = self.guess_residuetype(g[1])
            a = self.guess_elemtype(g[2])
            l = 'ATOM  {:5} {:3}  {:3}  {:4}   {:>8.3f} {:>8.3f} {:>8.3f}'.fromat(
                i+1,a,r[:3],g[1],g[4],g[5],g[6]
            )
            pdb += l + '\n'
        pdb += 'END\n'
        mol = Chem.MolFromPDBBlock(pdb,removeHs=False,sanitize=False)
        if 'moleculetype' in self.sections:
            name = self.sections['moleculetype'][0][0]
            mol.SetProp('_Name',name)
        mol.SetProp('Title',self.kws['title']) if 'title' in self.kws else mol.SetProp('Title','GMXTop')
        for i,a in enumerate(mol.GetAtoms()):
            a.SetProp('name',section[i][2])
        return mol

    def get_mol_from_xyz(self,section=None):
        if not section:
            if not 'xyz' in self.sections: return
            section = self.sections['xyz']
        xyz = '{:}\nGMXTop\n'.format(len(section))
        for g in section:
            a = self.guess_elemtype(g[0])
            l = '{:}  {:}  {:}  {:}'.format(a,g[1],g[2],g[3])
            xyz += l + '\n'
        xyz += '\n\n'
        mol = Chem.MolFromXYZBlock(xyz,removeHs=False,sanitize=False)
        if 'moleculetype' in self.sections:
            name = self.sections['moleculetype'][0][0]
            mol.SetProp('_Name',name)
        mol.SetProp('Title',self.kws['title']) if 'title' in self.kws else mol.SetProp('Title','GMXTop')
        for i,a in enumerate(mol.GetAtoms()):
            a.SetProp('name',section[i][0])
        return mol

    def guess_elemtype(self,s):
        s = s.strip()
        if not s:
            return '*'
        if s[:2] in self._periodic_table:
            return s[:2]
        if s[0] in self._periodic_table:
            return s[0]
        return s

    def guess_residuetype(self,s):
        s = s.strip()
        if not s:
            return 'UNK'
        if len(s) <= 2:
            return s
        aa = [
            'ASP', 'ALA', 'ARG', 'ASN', 'CYS', 'GLN', 'GLU', 'GLY', 'ILE', 'LEU',
            'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HIS',
        ]
        if s.upper() in aa:
            return s.upper()
        if s[:3].upper() in aa:
            return s[:3].upper()
        return s

    def file_to_pdb(self,file=None):
        mol = self.get_mol()
        if not mol:
            print('Warning: cannot generate PDB file')
            return
        new = self.genfnew(file,ext='.pdb')
        Chem.MolToPDBFile(mol,new)

    def file_to_gro(self,file=None):
        if 'gro' in self.sections:
            grosection = self.sections['gro']
        elif 'atoms' in self.sections:
            mol = self.get_mol_from_top()
            if not mol:
                print('Warning: cannot generate GRO file')
                return
            grosection = []
            conf = mol.GetConformer()
            g = self.sections['atoms']
            for i,p in enumerate(conf.GetPositions()):
                grosection.append((g[i][2],g[i][3],g[i][4],p[0],p[1],p[2]))
        elif 'xyz' in self.sections:
            grosection = []
            g = self.sections['atoms']
            for i,p in enumerate(self.sections['xyz']):
                grosection.append((g[i][2],g[i][3],g[i][4],p[0],p[1],p[2]))
        else:
            print('Warning: cannot generate GRO file')
            return

        if 'moleculetype' in self.sections:
            out = self.sections['moleculetype'][0][0]
        else:
            out = self.sections['title'] if 'title' in self.kws else 'GMXTopParser'
        out += '\n' + str(len(grosection)) + '\n'
        for g in grosection:
            out += '{:5}{:<5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}'.format(*g[:7])
            if len(g) > 7:
                out += '{:8.4f}{:8.4f}{:8.4f}'.format(g[7],g[8],g[9])
            out += '\n'
        if 'grobox' in self.kws:
            out += self.kws['grobox']
        else:
            xl = min([i[4] for i in grosection])
            xh = max([i[4] for i in grosection])
            yl = min([i[5] for i in grosection])
            yh = max([i[5] for i in grosection])
            zl = min([i[6] for i in grosection])
            zh = max([i[6] for i in grosection])
            l = '{:8.3f}{:8.3}{:8.3f}'.format(xh-xl,yh-yl,zh-zl)
            out += l + '\n'
            self.kws['grobox'] = l

        new = self.genfnew(file,ext='.gro')
        with open(new,'wt') as f:
            f.write(out)

    def file_to_itp(self,file=None):
        new = self.genfnew(file,ext='.itp')
        with open(new,'wt') as f:
            f.write(str(self))

    def genfnew(self,file=None,ext=None):
        if file:
            if os.path.isfile(file): return file
            file = '_'.join(file.split())
        else:
            file = 'gmx_base'
        if ext and not file.endswith(ext):
            file += ext
        i = 1
        while True:
            new = 'gmx-'+str(i)+'-'+file
            if not os.path.isfile(new): break
            i += 1
        print(f'Note: generating file: {new}')
        return new


if __name__ == '__main__':
    file = 'SystemFull.top'
    parser = GromacsTopParser()
    top = parser.read(file)
    op = GromacsTopOps(top=top)
    op.get_mol_from_top(2)




