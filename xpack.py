#!/usr/bin/env python3

import os
import sys
import io
import re
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
    'version 0.6.0 : add `pyminifier`',
    'version 0.7.0 : add `--mini-files`',
    'version 0.8.0 : add `--show-source-stdout`',
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


def xuntokenize(tokens):
    """tokens: [[tokentype, string, (srow, scol), (erow, ecol, line)], ...]"""
    out = ""
    last_lineno = -1
    last_col = 0
    for tok in tokens:
        token_string = tok[1]
        start_line, start_col = tok[2]
        end_line, end_col = tok[3]
        # The following two conditionals preserve indentation:
        if start_line > last_lineno:
            last_col = 0
        if start_col > last_col and token_string != '\n':
            out += (" " * (start_col - last_col))
        out += token_string
        last_col = end_col
        last_lineno = end_line
    return out


def xtokenize(source):
    """Tokenizes *source* and returns the tokens as a list of lists"""
    if isinstance(source,str):
        if os.path.isfile(source):
            with open(source,'rt') as f:
                source = f.read()
    else:
        print('Fatal: tokenizer: only file or string input supported')
        return []
    sio = io.StringIO(source)
    return [list(a) for a in tokenize.generate_tokens(sio.readline)]


def remove_blank_lines_and_spaces(source):
    bo = False
    if isinstance(source,(list,tuple)):
        source = xuntokenize(source)
        bo = True
    result = ''
    for line in source.split('\n'):
        l = line.rstrip()
        if len(l):
            result += l + '\n'
    if bo:
        return xtokenize(result)
    return result


def remove_blank_lines_and_spaces2(source):
    """to use, comments should be removed firstly"""
    bo = False
    if isinstance(source,(list,tuple)):
        source = xuntokenize(source)
        bo = True
    result = ''
    for line in source.split('\n'):
        result += line.rstrip() + '\n'              # only remove right whitespace
    tokens = xtokenize(result)
    prev = -1
    for tok in tokens:
        if tok[0] == 61:
            if prev == 4:   # newline already exist
                tok[1] = ''
                tok[3] = (tok[2][0], tok[2][1])      # tuple
        else:
            prev = tok[0]
    if bo:
        return tokens
    return xuntokenize(tokens)


def remove_comments(tokens):
    for tok in tokens:
        if tok[0] == tokenize.COMMENT:
            tok[1] = ''
            tok[3] = (tok[2][0], tok[2][1])      # tuple


def remove_docstrings(tokens):
    """must be docstring, for token type:

    > token.STRING: 3
    > token.NEWLINE: 4
    > token.INDENT: 5
    > token.DEDENT: 6
    > token.EXCLAMATION: 54
    > token.FSTRING_START: 61       # pure new line, changed from "NL"
    > token.NT_OFFSET: 256

    1) (filebeginning), {3}
    2) 61 -> {3}
    3) 6  -> {3}
    4) 54 -> 4 -> [61] -> 5 -> {3}  # `61` is optional, can be multiple
    """
    if not tokens: return tokens
    sall = ''.join(['{:3}'.format(i[0]) for i in tokens])
    sall = sall.replace(' ','@')
    if sall[:3] == '@@3':
        tokens[0][1] = ''
        tokens[0][3] = (tokens[0][3][0],  tokens[0][2][1])      # tuple
    # search pattern: @@4@@3, @61@@3, @@6@@3, @54@@4{@61}*@@5@@3
    p = re.compile(r'@@4@@3|@61@@3|@@6@@3|@54@@4(@61)*@@5@@3')
    for i in p.finditer(sall):
        k = i.span()[1] // 3 - 1        # starts from 1
        tokens[k][1] = ''
        tokens[k][3] = (tokens[k][2][0],  tokens[k][2][1])      # tuple


def fix_empty_methods(tokens):
    """
    empty method will become:
        def func(): pass
        class cls(): pass
    """
    bo = False
    if isinstance(tokens,(list,tuple)):
        tokens = xuntokenize(tokens)
        bo = True
    lines = tokens.split('\n')
    if not lines[-1]:
        lines.pop(len(lines)-1)
    if not lines: return ''
    for k in range(len(lines)-1):
        if not lines[k]: continue       # inline docstring may have empty line
        if lines[k][-1] == ':':
            p = lines[k].lstrip()       # may use `\t` or spaces
            if p[:3] == 'def' or p[:5] == 'class':
                offset = len(lines[k]) - len(p)
                u = lines[k+1].lstrip()
                offnew = len(lines[k+1]) - len(u)
                if offset == offnew and (u[:3] == 'def' or u[:5] == 'class'):
                    lines[k] += ' pass'
    p = lines[-1].lstrip()
    if (p[:3] == 'def' or p[:5] == 'class') and p[-1] == ':':
        lines[-1] += ' pass'
    result = '\n'.join(lines)
    if bo:
        return xtokenize(result)
    return result


def reduce_operators(source):
    """
    Remove spaces between operators in *source* and returns the result.
    Example::

        def foo(foo, bar, blah):
            test = "This is a %s" % foo

    Will become::

        def foo(foo,bar,blah):
            test="This is a %s"%foo

    Trailing commas cannot be removed, cause we do not know how is used:
    e.g.;
        for i in (1,): print(i)         # this makes int 1 iterable
    """
    io_obj = io.StringIO(source)
    prev_tok = [-99,]            # xzdebug, for the beginning of file `__doc__`
    out = ""
    last_lineno = -1
    last_col = 0
    nl_types = (tokenize.NL, tokenize.NEWLINE)
    joining_strings = False
    new_string = ""
    for tok in tokenize.generate_tokens(io_obj.readline):
        token_type = tok[0]
        token_string = tok[1]
        start_line, start_col = tok[2]
        end_line, end_col = tok[3]
        if start_line > last_lineno:
            last_col = 0
        if token_type != tokenize.OP:
            if start_col > last_col and token_type not in nl_types:
                if prev_tok[0] != tokenize.OP:
                    out += (" " * (start_col - last_col))
            if token_type == tokenize.STRING:
                if prev_tok[0] == tokenize.STRING:
                    # Join the strings into one
                    string_type = token_string[0] # '' or ""
                    prev_string_type = prev_tok[1][0]
                    out = out.rstrip(" ") # Remove any spaces we inserted prev
                    if not joining_strings:
                        # Remove prev token and start the new combined string
                        out = out[:(len(out)-len(prev_tok[1]))]
                        prev_string = prev_tok[1].strip(prev_string_type)
                        new_string = (
                            prev_string + token_string.strip(string_type))
                        joining_strings = True
                    else:
                        new_string += token_string.strip(string_type)
        else:
            # if token_string in ('}', ')', ']'):           # xzdebug
            #     if prev_tok[1] == ',':
            #         out = out.rstrip(',')
            if joining_strings:
                # NOTE: Using triple quotes so that this logic works with
                # mixed strings using both single quotes and double quotes.
                out += "'''" + new_string + "'''"
                joining_strings = False
            if token_string == '@': # Decorators need special handling
                if prev_tok[0] == tokenize.NEWLINE:
                    # Ensure it gets indented properly
                    out += (" " * (start_col - last_col))
        if not joining_strings:
            out += token_string
        last_col = end_col
        last_lineno = end_line
        prev_tok = tok
    return out


def join_multiline_pairs(source, pair="()"):
    """
    Idea is to remove `\n` in every defined pair `OP`

    Finds and removes newlines in multiline matching pairs of characters in
    *source*.

    By default it joins parens () but it will join any two characters given via
    the *pair* variable.

    .. note::

        Doesn't remove extraneous whitespace that ends up between the pair.
        Use `reduce_operators()` for that.

    Example::

        test = (
            "This is inside a multi-line pair of parentheses"
        )

    Will become::

        test = (            "This is inside a multi-line pair of parentheses"        )

    """
    opener = pair[0]
    closer = pair[1]
    out_tokens = []
    open_count = 0
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        token_type = tok[0]
        token_string = tok[1]
        if token_type == tokenize.OP and token_string in pair:
            if token_string == opener:
                open_count += 1
            elif token_string == closer:
                open_count -= 1
            out_tokens.append(tok)
        elif token_type in (tokenize.NL, tokenize.NEWLINE):
            if open_count == 0:
                out_tokens.append(tok)
        else:
            out_tokens.append(tok)
    return xuntokenize(out_tokens)


def dedent(source, use_tabs=False):
    """
    Minimizes indentation to save precious bytes.  Optionally, *use_tabs*
    may be specified if you want to use tabulators (\t) instead of spaces.

    Example::

        def foo(bar):
            test = "This is a test"

    Will become::

        def foo(bar):
         test = "This is a test"
    """
    indent_char = '\t' if use_tabs else ' '
    out = ""
    last_lineno = -1
    last_col = 0
    prev_start_line = 0
    indentation_level = 0
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        if tok[0] == tokenize.INDENT:
            indentation_level += 1
            continue
        if tok[0] == tokenize.DEDENT:
            indentation_level -= 1
            continue
        token_string = tok[1]
        start_line, start_col = tok[2]
        end_line, end_col = tok[3]
        if start_line > last_lineno:
            last_col = 0
        if start_line > prev_start_line:
            if token_string in (',', '.'):
                out += token_string
            else:
                indentation = indent_char * indentation_level
                out += indentation + token_string
        elif start_col > last_col:
            out += indent_char + token_string
        else:
            out += token_string
        prev_start_line = start_line
        last_col = end_col
        last_lineno = end_line
    return out


def minify(tokens, use_tabs=None):
    """tokens: file-contents or tokens"""
    if not isinstance(tokens, (list,tuple)):
        tokens = xtokenize(tokens)
    remove_comments(tokens)
    result = xuntokenize(tokens)
    result = remove_blank_lines_and_spaces2(result)
    result = join_multiline_pairs(result)
    result = join_multiline_pairs(result, '[]')
    result = join_multiline_pairs(result, '{}')
    tokens = xtokenize(result)
    remove_docstrings(tokens)
    result = xuntokenize(tokens)
    result = remove_blank_lines_and_spaces2(result)
    result = result.lstrip('\n')
    result = fix_empty_methods(result)
    result = reduce_operators(result)
    result = dedent(result, use_tabs)
    return result


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
            if isinstance(t[1],str):
                if (os.path.isfile(t[0]) and hasattr(t[0],'read')):
                    finals.append((t[0].read(),t[1]))
                elif isinstance(t[0],str):      # str, source contents
                    finals.append(t)
                else:
                    bo = True
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
        print('Note: use python source: zipfile::__main__')
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
            w.write(t[0])
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
    print(f'Note: overall size reduction: {p} of original size\n')


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

# check going to be processed files
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
        '-m', '--mini',
        action='store_true',
        help='minification by removing docstrings, comments, white spaces, empty lines, and indentions'
    )
    parser.add_argument(
        '-O', '--show-source-stdout',
        action='store_true',
        help='show sources to stdout, valid when `--mini` is used'
    )
    parser.add_argument(
        '-mf', '--mini-files',
        action='store_true',
        help='minification files only, but not create executable'
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
    parser.add_argument(
        '--features',
        action='store_true',
        help='show development features'
    )
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    if args.features:
        for i in FEATURES: print(i)
        sys.exit(0)

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

    if args.list_modules:
        print('\nNote: going to include files:')
        for i in files:
            print(f' ->  {i}')
        if args.mini:
            print('Note: minification will be performed')
        print()
        return

    if args.list_module_variables:
        for f in files:
            results = get_module_variables(f)
            print(f'\n>>>:: module variables for file: {f}')
            if results:
                print('  '.join(results))
        print()
        return

    if args.mini or args.mini_files:
        new = []
        for f in files:
            if f.endswith('.py'):
                n = minify(f)
                new.append((n,f))
            else:
                new.append((f,f))
        files = new
        print('\nNote: minification is performed')
    
        if args.show_source_stdout:
            for c,f in files:
                print('\n\n# ' + f)
                print(c,end='')
                print()
            print()
            return

    if args.mini_files:
        dest = '__xpackmini'
        os.makedirs(dest,exist_ok=True)
        pures = [t for t in files if not os.path.isfile(t[0])]
        cwd = os.getcwd()
        os.chdir(dest)
        for t in pures:
            base = os.path.basename(t[1])
            with open(base,'wt') as f: f.write(t[0])
        os.chdir(cwd)
        print(f'Note: check folder for results: {dest}\n')
        return

    bo = True if args.precompile else False
    _zip_packall(files,precompile=bo)


if __name__ == '__main__':
    main()




