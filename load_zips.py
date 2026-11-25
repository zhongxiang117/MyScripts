import io
import os
import sys
import zipfile
import importlib.abc
import importlib.util


class MergedZipSourceLoader(importlib.abc.SourceLoader):
    """
    Inheritance:
        importlib.abc.Loader
            ↓ inherits
        importlib.abc.InspectLoader         get_source | get_code
            ↓ inherits
        importlib.abc.SourceLoader          get_filename
            ↓ inherits
        importlib.abc.ResourceLoader        get_data
    """
    def __init__(self, origin, srcbytes):
        self.origin = origin
        self.srcbytes = srcbytes

    def get_filename(self, fullname):
        """__file__"""
        return self.origin

    def get_source(self, fullname):
        """contents"""
        return self.srcbytes.decode('utf-8')

    def get_code(self, fullname):
        """compiled object"""
        src = self.get_source(fullname)
        return compile(src, self.origin, 'exec')

    def get_data(self, path):
        """binary"""
        return self.srcbytes


class MergedZipFinder(importlib.abc.MetaPathFinder):
    """
    Assume:
        all.zip
        ├─ a.py
        ├─ dira/
        │  ├─ __init__.py
        │  ├─ b.py
        │  └─ diraa/
        │     ├─ __init__.py
        │     └─ bb.py
        └─ sub.zip
           ├─ s.py
           ├─ a.py
           ├─ pot.zip
           │  └─ c.py
           └─ dirb/
               ├─ __init__.py
               └─ b.py

    1) recursive_level == 0:
        only import codes in `all.zip`
    2) recursive_level == 1:
        import codes from `all.zip -> sub.zip`
    3) recursive_level == n:
        import `all.zip<0> -> sub.zip<1> -> pot.zip<2> -> ... -> *.zip<n>`

    Attributes:
        self.imports = [(entry, origin), ...]
        self.conflicts = {entry: set(origin), ...}
        self.backups = {origin: srcbytes, ...}
        self.recursive = {int: (subzipname, subzipbytes), ...}

        entry: means the content compressed inside zip, like `a.py`, folder `dira/`, or `sub.zip`
        origin: `<x:where_it_comes_from>`, like `<x:all.zip!a.py`, `<x:all.zip!sub.zip!s.py>`

    Usage:
        >>> finder = MergedZipFinder(zipfiles=['x.zip', 'y.zip'], recursive_level=2)
        >>> finder.run_recursively()
        >>> #finder.add_zipfile_bytes('w.zip')
        >>> sys.meta_path.insert(0, finder)

        >>> finder = MergedZipFinder()
        >>> finder.add_zipfile_bytes('x.zip')
        >>> find.load_from_bytes(filebytes, 'nice_module')
        >>> sys.meta_path.insert(0, finder)

    Call Sequence:
        sys.meta_path -> find_spec
                            => spec_from_loader => loader.exec_module ~> get_code ~> get_source
    """
    def __init__(self, zipfiles=None, recursive_level=0):
        self.imports = []
        self.conflicts = {}
        self.backups = {}
        self.recursive = {}
        self.zipfiles = []
        if isinstance(zipfiles, (str, bytes)): zipfiles = [zipfiles, ]
        if isinstance(zipfiles, (tuple, list)):
            for f in zipfiles:
                if isinstance(f, bytes) or os.path.isfile(f):
                    self.zipfiles.append(f)
        self.recursive_level = recursive_level

    def run_recursively(self):
        for f in self.zipfiles:
            news = self.add_zipfile_bytes(f)
            self.recursive.setdefault(0, []).extend(news)
        for i in range(0, self.recursive_level):
            if i not in self.recursive: continue
            news = self.recursive[i]
            for g in news:
                sub = self.add_zipfile_bytes(g[1], label=g[0])
                self.recursive.setdefault(i+1, []).extend(sub)

    def add_zipfile_bytes(self, zipbytes, label=None):
        newzipbytes = []
        if os.path.isfile(zipbytes):
            bio = zipbytes
            if not label: label = zipbytes
        else:
            bio = io.BytesIO(zipbytes)
            if not label: label = 'xbytes'
        with zipfile.ZipFile(bio, 'r') as zf:
            for o in zf.infolist():
                entry = o.filename
                if entry.endswith('/'):
                    continue

                if entry.endswith('.zip'):
                    newlabel = '{:}!{:}'.format(label,entry)
                    newbytes = zf.read(entry)
                    newzipbytes.append((newlabel, newbytes))
                    continue

                if not entry.endswith('.py'): continue

                srcbytes = zf.read(entry)
                origin = '<x:{:}!{:}>'.format(label, entry)
                self.backups[origin] = srcbytes

                use = [s[1] for s in self.imports if s[0] == entry]
                if use:
                    self.conflicts.setdefault(entry, set()).update((origin, *use))
                    c = ' | '.join(use)
                    print('Conflict: {:}: exist {:}'.format(origin, c))
                    continue

                self.imports.append((entry,origin))

        return newzipbytes

    def find_spec(self, fullname, path, target=None):
        base = fullname.replace('.', '/')
        candidates = [(f'{base}.py', False), (f'{base}/__init__.py', True)]
        for entry, is_pkg in candidates:
            use = [s[1] for s in self.imports if s[0] == entry]
            if not use: continue
            origin = use[0]
            srcbytes = self.backups[origin]
            loader = MergedZipSourceLoader(origin, srcbytes)
            spec = importlib.util.spec_from_loader(fullname, loader, origin=origin, is_package=is_pkg)
            if is_pkg:
                spec.submodule_search_locations = [origin]
            return spec
        return None

    def load_from_bytes(self, srcbytes, module, force=False):
        """load srcbytes into `sys.modules` as the identifier `module`,
            whether force overwrite when exists
        """
        if module in sys.modules and not force:
            print('Warning: already exist: {:}'.format(module))
            return
        loader = MergedZipSourceLoader(module, srcbytes)
        spec = importlib.util.spec_from_loader(module, loader, origin=module, is_package=False)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        sys.modules[module] = mod



