#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    tenxtools.svgfigure
    ~~~~~~~~~~~~~~~~~~~

    SVG Figure

    @Copyright: (c) 2018-09 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM
import svgwrite

from bioutensil.bplot import svgbase
from bioutensil.bplot import text


class SVGFigure():

    def __init__(self,
                 out_prefix='demo',
                 format='svg',
                 title='SVG Figure',
                 size=(1000, 1000),
                 version='0.0', author='Lingxi Chen', email='chanlingxi@gmail.com', address="City University of Hong Kong",
                 info={'Info': 'Please tell something about this figure.'},
                 icon_related_insert='bottom:middle', info_related_insert='up:middle',
                 **kwargs
                 ):
        self.out_prefix = out_prefix
        self.format = format
        self.title = title
        self.size = size

        self.version = version
        self.author = author
        self.email = email
        self.address = address
        self.info = info
        self.info_related_insert = info_related_insert
        self.icon_related_insert = icon_related_insert

        self.dwg = svgwrite.Drawing(out_prefix + '.svg', size=size)
        self.canvas = svgbase.Container((0, 0), size).svgobj

    def draw(self):

        self._draw_info()
        self._draw_icon()

        self.dwg.add(self.canvas)
        self.dwg.save()
        self._convert()

    def _convert(self):
        if self.format == 'pdf':
            drawing = svg2rlg(self.out_prefix + '.svg')
            renderPDF.drawToFile(drawing, '{}.pdf'.format(self.out_prefix))

        if self.format == 'png':
            drawing = svg2rlg(self.out_prefix + '.svg')
            renderPM.drawToFile(drawing, "{}.png".format(self.out_prefix))

    def _draw_textarea(self, text_str, insert, parent, related_insert):
        textarea = text.TextArea(text_str, insert,
                                 parent=parent, related_insert=related_insert).render()
        return textarea

    def _draw_icon(self):
        icon_text = '\n'
        icon_text += '{}\n'.format(self.title)
        icon_text += 'v{}\n'.format(self.version)
        icon_text += 'Author: {}\n'.format(self.author)
        icon_text += 'Email: {}\n'.format(self.email)
        icon_text += self.address
        icon_text += '\n'
        icon = self._draw_textarea(icon_text, (0, 0),
                                   parent=self.canvas, related_insert=self.icon_related_insert)
        self.canvas.add(icon)

    def _draw_info(self):
        infos = ['', self.title]
        infos += [name + ': ' + value for name, value in self.info.items()]
        infos += ['']
        info_text = '\n'.join(infos)
        info = self._draw_textarea(info_text, (0, 0),
                                   parent=self.canvas, related_insert=self.info_related_insert)
        self.canvas.add(info)
