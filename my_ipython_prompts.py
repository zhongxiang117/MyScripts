from IPython.terminal.prompts import Prompts, Token
from platform import python_version

import os
import sys

VERSION = python_version()
CONDA = os.getenv('CONDA_DEFAULT_ENV','')
LEVEL = os.getenv('SHLVL','')

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


ip = get_ipython()
ip.prompts = MyPrompt(ip)

#if 'schrodinger' in sys.prefix.lower():
#    title = 'level-%s: Python %s: Schrodinger: {cwd}' % (LEVEL, VERSION)
#else:
#    title = 'level-%s: Python %s: {cwd}' % (LEVEL, VERSION)
#ip.prompts.shell.term_title_format = title

if 'schrodinger' in sys.prefix.lower():
    ip.prompts.shell.term_title_format = 'Schrodinger: {cwd}'



