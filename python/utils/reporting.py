#!/usr/bin/env python
"""
Module holding some definitions for automatically creating latex reports
"""

from contextlib import contextmanager
from os import environ


BASIC_PREAMBLE = r'''\documentclass[a4paper,11pt]{scrartcl}
\usepackage{graphicx}
\usepackage{subfig}
\usepackage{amsmath}
\usepackage[margin=2cm]{geometry}
\usepackage{tabulary}

\newcommand{\rootfig}[2]{\includegraphics[width=#1\linewidth, angle=270, origin=c]{#2}}

'''

try:
    LATEX_EXE = environ['LATEX_EXE']
except KeyError:
    LATEX_EXE = 'pdflatex'

PLOT_COMMAND = r'\includegraphics[width=0.5\linewidth]'
if LATEX_EXE == 'xelatex':
    PLOT_COMMAND = r'\rootfig{0.475}'

MAX_FIG_P_PAGE = 6

@contextmanager
def open_tex_file(filename, preamble=BASIC_PREAMBLE):
    """
    Context manager to open a tex file and populate it with a (customizable
    preamble). Will also add the begin/end of document statements before/after
    the content written to the file.

    Args:
        filename (str): The output filename of the tex file
        preamble (str, optional): The preamble that should be used for the tex
            file

    Yields:
        file object: File object that can be written to. Everything written to
            this file will be considered content.
    """
    with open(filename, 'w') as texfile:
        texfile.write(preamble)
        texfile.write(r'\begin{document}')
        texfile.write('\n\n')

        yield texfile

        texfile.write('\n\n')
        texfile.write(r'\end{document}')


def create_subfloat(plot, label):
    """
    Create a subfloat entry
    """
    return (r'\subfloat[][LABEL]{PLTCMD{PLOT}}'
            .replace('LABEL', label).replace('PLOT', plot)
            .replace('PLTCMD', PLOT_COMMAND))


def create_figure(plots):
    """
    Make the latex figure
    """
    ret_str = []

    fig_started = False
    for iplot, (label, pname) in enumerate(plots.items()):
        if iplot % MAX_FIG_P_PAGE == 0:
            if fig_started:
                ret_str.append(r'\end{figure}')
                ret_str.append('')
                fig_started = False
            if not fig_started:
                ret_str.append(r'\begin{figure}')
                ret_str.append('')
                fig_started = True

        ret_str.append(create_subfloat(pname, label))

        if iplot % 2 != 0:
            ret_str.append('')

    if fig_started:
        ret_str.append(r'\end{figure}')

    return '\n'.join(ret_str)
