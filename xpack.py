#!/usr/bin/env python3

import os
import sys
import argparse


__version__ = '0.1.0'
__author__ = 'Xiang Zhong'


def enumerate_local_modules(path,include_sofile=None,strip=False):
    """return all python source and soname in `path`"""
    local_modules = set()
    old = os.getcwd()
    full = os.path.abspath(path)
    cwd = os.path.basename(full) if os.path.isfile(full) else full
    os.chdir(cwd)
    # Now check the local dir for matching modules
    for root, dirs, files in os.walk('./'):
        for f in files:
            bo = False
            if f.endswith('.py'):
                bo = True
            if include_sofile and f.endswith('.so'):
                bo = True
            if bo:
                if strip:
                    f = f[:-3]
                module = os.path.join(root,f)
                local_modules.add(module)
    os.chdir(old)
    return list(local_modules)


def _zip_packall(fileobjs):
    """
    Args:
        fileobjs : [(fobj, modulename), filename, ...]      # tuple 2 | file

    * The first item will be used as `__main__.py`
    """
    import zipfile
    import tempfile
    finals = []
    for t in fileobjs:
        bo = False
        if isinstance(t,str) and os.path.isfile(t):
            finals.append((t,t))
        elif isinstance(t,(list,tuple)) and len(t) == 2:
            if (os.path.isfile(t[0]) or hasattr(t[0],'read')) and isinstance(t[1],str):
                finals.append(t)
            else:
                bo = True
        else:
            bo = True
        if bo:
            print(f'Warning: zip_py: not valid: ignoring: {t}')
    if not finals:
        print('Fatal: zip_py: no inputs')
        return

    # This is so it will still execute as a zip
    finals[0] = (finals[0][0],'__main__.py')

    i = 1
    while True:
        dest = 'xz' + str(i) + '.zip'
        if not os.path.isfile(dest):
            break
        i += 1
    print(f'Note: zipfile will be save to: {dest}')
    # Hopefully some day we'll be able to use ZIP_LZMA to save even more space
    z = zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED)
    total_old_size = 0.0
    for t in finals:
        if isinstance(t[0],str):
            file = t[0]
        else:
            w = tempfile.NamedTemporaryFile(mode='w')
            w.write(t[0].read())
            w.flush()
            file = w.name
        total_old_size += os.path.getsize(file)
        print(f'Note: zipfile writing::: {t[1]}')
        z.write(file, t[1])
        if not isinstance(t[0],str):
            w.close()
    z.close()

    i = 1
    while True:
        exe = 'xp' + str(i)
        if not os.path.isfile(exe):
            break
        i += 1
    print(f'Note: executable will be save to: {exe}')

    with open(exe, 'wb') as f:
        f.write(b'#!/usr/bin/env python3\n')
        with open(dest, 'rb') as r:
            f.write(r.read())
    # Make it executable since shebang added
    os.chmod(exe, 0o755)
    p = round(os.path.getsize(exe)/total_old_size, 4)
    print(f'Note: overall size reduction: {p} of original size')


def main():
    parser = argparse.ArgumentParser(
        usage="""xpyminifier:

# pack all files
>>> xpack.py -s scrfile1 [scrfile2, ...]  [--recursive]

""",
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=__version__
    )
    parser.add_argument(
        '-s', '--srcfiles',
        dest='srcfiles',
        nargs='+',
        help='source files need to be processed'
    )
    parser.add_argument(
        '-R', '--recursive',
        dest='recursive',
        action='store_true',
        help='input `-s` is a folder, recursively search all its modules'
    )
    parser.add_argument(
        '-S', '--include_sofile',
        dest='include_sofile',
        action='store_true',
        help='include sona files, valid when `-R` used'
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    if args.srcfiles:
        if not os.path.isfile(args.srcfiles[0]):
            print(f'Fatal: the first input is not a file: {args.srcfiles[0]}')
    else:
        print('Fatal: no inputs')
        sys.exit(0)

    files = [i for i in args.srcfiles if os.path.isfile(i)]
    if args.recursive:
        for i in args.srcfiles:
            if os.path.isdir(i):
                m = enumerate_local_modules(i,include_sofile=args.include_sofile,strip=False)
                for j in m:
                    files.append(os.path.join(i,j))

    _zip_packall(files)


if __name__ == '__main__':
    main()




