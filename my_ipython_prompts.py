from IPython.terminal.prompts import Prompts, Token
#from platform import python_version

import os
import sys

FEATURES = [
    'version 0.1.0  : for IPython Prompts',
    'version 0.2.0  : add parser for ``help``',
    'version 0.3.0  : more versatile on ``help``',
    'version 0.4.0  : make more robust on ``help``',
]

__version__ = FEATURES[-1].split(':')[0].split()[1]

#VERSION = python_version()
CONDA = os.getenv('CONDA_DEFAULT_ENV','')
#LEVEL = os.getenv('SHLVL','')
#if 'schrodinger' in sys.prefix.lower():
#    title = 'level-%s: Python %s: Schrodinger: {cwd}' % (LEVEL, VERSION)
#else:
#    title = 'level-%s: Python %s: {cwd}' % (LEVEL, VERSION)


class MyPrompt(Prompts):
    def in_prompt_tokens(self):
        return [
            (Token.Prompt, self.vi_mode() ),
            (Token.Prompt, '('+CONDA+') '),
            (Token.Literal.String, os.path.basename(os.getcwd())+': '),
            (Token.Prompt, 'In ['),
            (Token.PromptNum, str(self.shell.execution_count)),
            (Token.Prompt, ']: '),
        ]


def my_cmd_help(lines):
    news = []
    for line in lines:
        if not (line.startswith('%') or line.startswith('help(')):
            tmp = line.rstrip()
            if tmp.endswith(' %h'):
                tmp = tmp[:-2]
                bo = True
            elif tmp.endswith(' %help'):
                tmp = tmp[:-5]
                bo = True
            else:
                bo = False
            if bo:
                # cases (simplified):
                ##  os.getenv %h
                ##  os.getenv('PATH') %h
                ##  v = os.getenv('PATH') %h
                ndx = tmp.find('(')
                if ndx != -1:
                    tmp = tmp[:ndx]
                ndx = tmp.find('=')
                if ndx != -1:
                    tmp = tmp[ndx+1:]
                line = 'help('+tmp+')'
        news.append(line)
    return news


ip = get_ipython()
ip.prompts = MyPrompt(ip)

#ip.prompts.shell.term_title_format = title
if 'schrodinger' in sys.prefix.lower():
    ip.prompts.shell.term_title_format = 'Schrodinger: {cwd}'


ip.input_transformers_cleanup.append(my_cmd_help)




