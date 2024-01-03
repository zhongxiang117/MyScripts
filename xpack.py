#!/usr/bin/env python3

import os
import sys
import argparse
import tokenize
import py_compile


FEATURES = [
    'version 0.1.0 : start',
    'version 0.2.0 : add local modules filter for `builtins`',
    'version 0.3.0 : add option `--no-precompile` for `pyc` support',
    'version 0.4.0 : change option `--no-precompile` to `precompile`',
    'version 0.5.0 : add parser for non-py file',
    'version 0.5.1 : avoid duplicate file',
]

__version__ = FEATURES[-1].split()[1]
__author__ = 'Xiang Zhong'

builtins = [
    'abs',  'all',  'any',  'bin',  'bool',  'bytearray',  'bytes',  'callable',
    'chr',  'classmethod',  'compile',  'complex',  'delattr',  'dict',  'dir',
    'divmod',  'enumerate',  'eval',  'filter',  'float',  'format',  'frozenset',
    'getattr',  'globals',  'hasattr',  'hash',  'hex',  'id',  'input',  'int',
    'isinstance',  'issubclass',  'iter',  'len',  'list',  'locals',  'map',
    'max',  'min',  'next',  'object',  'oct',  'open',  'ord',  'pow',  'print',
    'property',  'quit',  'range',  'raw_input',  'repr',  'reversed',  'round',
    'set',  'setattr',  'slice',  'sorted',  'staticmethod',  'str',  'sum',
    'super',  'tuple',  'type',  'vars',  'zip',

    'os', 'sys', 'argparse', 'import',
    'is', 'None', 'if', 'else', 'elif', 'and', 'or',
]

def enumerate_local_modules(path,include_files=None,strip=False):
    """return all python source and include_files in `path`"""
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
            if include_files and f in include_files:
                bo = True
            if bo:
                if strip:
                    f = f[:-3]
                module = os.path.join(root,f)
                local_modules.add(module)
    os.chdir(old)
    return list(local_modules)


def _zip_packall(fileobjs,precompile=True):
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
    if finals[0][0].endswith('.py'):
        finals[0] = (finals[0][0],'__main__.py')
    elif finals[0][0].endswith('.pyc'):
        finals[0] = (finals[0][0],'__main__.pyc')
    else:
        print('Warning: not a python source: zipfile::__main__')
        finals[0] = (finals[0][0],'__main__.py')

    i = 1
    while True:
        dest = 'xz' + str(i) + '.zip'
        if not os.path.isfile(dest):
            break
        i += 1
    print(f'Note: zipfile will be saved to: {dest}')
    # Hopefully some day we'll be able to use ZIP_LZMA to save even more space
    z = zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED)
    total_old_size = 0.0
    for t in finals:
        if os.path.isfile(t[0]):
            file = t[0]
        else:
            w = tempfile.NamedTemporaryFile(mode='w',suffix='.py')
            w.write(t[0].read())
            w.flush()
            file = w.name
        total_old_size += os.path.getsize(file)
        if precompile and file.endswith('.py'):
            file = py_compile.compile(file,file+'c')
            m = t[1] if t[1].endswith('.pyc') else t[1]+'c'
        else:
            m = t[1]
        print(f'Note: zipfile writing::: {m}')
        z.write(file, m)
        if not os.path.isfile(t[0]):
            w.close()
    z.close()

    i = 1
    while True:
        exe = 'xp' + str(i)
        if not os.path.isfile(exe):
            break
        i += 1
    print(f'Note: executable will be saved to: {exe}')
    s = 'Has' if precompile else 'No'
    print(f' ->  {s} precompile used')

    with open(exe, 'wb') as f:
        f.write(b'#!/usr/bin/env python3\n')
        with open(dest, 'rb') as r:
            f.write(r.read())
    # Make it executable since shebang added
    os.chmod(exe, 0o755)
    p = round(os.path.getsize(exe)/total_old_size, 4)
    print(f'Note: overall size reduction: {p} of original size')


def get_module_variables(file):
    if not file or not (isinstance(file,str) and os.path.isfile(file) and file.endswith('.py')):
        return []
    results = set()
    indent = 0
    operator = False
    prevstr = ''
    for tok in tokenize.generate_tokens(open(file,'rt').readline):
        if tok[0] == tokenize.INDENT:
            indent += 1
        elif tok[0] == tokenize.DEDENT:
            indent -= 1
        elif tok[0] == tokenize.OP:
            if tok[1] in ['(', '[','{']:
                operator = True
            elif tok[1] in [')', ']', '}']:
                operator = False
        elif tok[0] == tokenize.NAME and indent == 0 and not operator and prevstr != '.':
            results.add(tok[1])
        prevstr = tok[1]

    both = results.intersection(builtins)
    diff = results.difference(builtins)
    return sorted(list(diff))+sorted(list(both))


def main():
    parser = argparse.ArgumentParser(
        usage="""xpack.py:

# check files will be processed
>>> xpack.py -s scrfile1 [scrfile2, dir, ...]  -i add1.dat add2.md -L

# pack all files
>>> xpack.py -s scrfile1 [scrfile2, dir, ...]  -i add1.dat add2.md

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
        help='source files or folders need to be processed'
    )
    parser.add_argument(
        '-C', '--precompile',
        dest='precompile',
        action='store_true',
        help='precompile python source before writing to zip executable'
    )
    parser.add_argument(
        '-i', '--include-files',
        dest='include_files',
        nargs='+',
        help='include additional files, such as README or data file'
    )
    parser.add_argument(
        '-L', '--list-modules',
        dest='list_modules',
        action='store_true',
        help='tool, only list python module type variables'
    )
    parser.add_argument(
        '-M', '--list-module-variables',
        dest='list_module_variables',
        action='store_true',
        help='tool, only list python module type variables'
    )
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    if args.srcfiles:
        if not os.path.isfile(args.srcfiles[0]):
            print(f'Fatal: the first input is not a file: {args.srcfiles[0]}')
            sys.exit(0)
    else:
        print('Fatal: no inputs')
        sys.exit(0)

    files = []
    for i in args.srcfiles:
        if os.path.isfile(i):
            v = os.path.normpath(i)
            if v not in files:
                files.append(v)
    for i in args.srcfiles:
        if os.path.isdir(i):
            m = enumerate_local_modules(i,include_files=args.include_files,strip=False)
            for j in m:
                k = os.path.normpath(os.path.join(i,j))
                if k not in files:
                    files.append(k)
    bo = True
    if args.list_modules:
        print('\nNote: going to include files:')
        for i in files:
            print(f' ->  {i}')
        print()
        bo = False
    if args.list_module_variables:
        for f in files:
            results = get_module_variables(f)
            if results:
                print(f'\n>>>:: module variables for file: {f}')
                print('  '.join(results))
        print()
        bo = False
    if bo:
        bo = True if args.precompile else False
        _zip_packall(files,precompile=bo)


if __name__ == '__main__':
    main()




