#!/usr/bin/env python
"""
Module holding some definitions for automatically creating latex reports
"""

from contextlib import contextmanager


BASIC_PREAMBLE = r'''\documentclass[a4paper,11pt]{scrartcl}
\usepackage{graphicx}
\usepackage{subfig}
\usepackage{amsmath}
\usepackage[margin=2cm]{geometry}

\newcommand{\rootfig}[2]{\includegraphics[width=#1\linewidth, angle=270, origin=c]{#2}}

'''

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
